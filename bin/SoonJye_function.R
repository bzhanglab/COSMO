library(missForest)
library(biomaRt)
library(glmnet)
library(caret)
library(parallel)
library(doParallel)

run_2b <- function(pro_file, rna_file, anno_file, gene_file, out_dir="./",
                   clinical_attributes=NA){
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ########### Data Retrieval
  cat('Correcting dataset:', pro_file, " and ", rna_file , '... \n')
  out_files <- preprocess(rna_file, pro_file, gene_file, out_dir = out_dir)
  
  clinic  <- getClinical(anno_file)
  
  rnaseq <- getRNAseq(out_files[1], out_files[3])
  rna_atsm <- rnaseq$rna_atsm
  rna_sex  <- rnaseq$rna_sex
  
  proteome <- getProteome(out_files[2], out_files[3])
  pro_atsm <- proteome$pro_atsm
  pro_sex  <- proteome$pro_sex
  
  sample_n <- nrow(rna_atsm)
  sample_labels <- rownames(rna_atsm)
  cat('\nTotal', sample_n, 'samples... \n\n')
  rownames(rna_atsm) <- paste0('RNA_', 1:sample_n)
  rownames(rna_sex) <- paste0('RNA_', 1:sample_n)
  rownames(pro_atsm) <- paste0('PRO_', 1:sample_n)
  rownames(pro_sex) <- paste0('PRO_', 1:sample_n)
  
  
  
  ########### Screening RNA vs PRO mismatched samples with correlation
  #### first round of getting correlated genes
  cat('Screening First Round: ')
  corsample <- computeCorr(rna_atsm, pro_atsm)
  rankdist <- computeRankdist(corsample)
  
  pairdist <- getMatching(rankdist)
  nonmatch <- which(1:sample_n != pairdist$pro)
  cat('  Found', length(nonmatch), 'mismatched samples \n')
  
  #### second round of getting correlated genes (after removing nonmatch from first round)
  cat('Screening Second Round: ')
  corsample <- computeCorr(rna_atsm, pro_atsm, nonmatch=nonmatch)
  rankdist <- computeRankdist(corsample)
  pairdist <- getMatching(rankdist)
  
  rnamatch <- pairdist$pro
  promatch <- pairdist$rna[order(pairdist$pro)]
  nonmatch <- which(1:sample_n != rnamatch)
  cat('  Found', length(nonmatch), 'mismatched samples \n')
  
  
  rownames(corsample) <- sample_labels
  colnames(corsample) <- sample_labels
  rownames(rankdist)  <- sample_labels
  colnames(rankdist)  <- sample_labels
  
  # output correlation file
  cor_file <- paste0(out_dir,"/correlation.csv")
  write.table(corsample, cor_file , col.names=TRUE, row.names=TRUE, sep=',')
  
  #### Determine spurious match due to being left out (threshold is log of n)
  spumatch <- which(1:sample_n == rnamatch & pairdist$distance >= log(sample_n))
  if (length(spumatch) > 0){
    nonmatch <- c(nonmatch, spumatch)
    cat('  Spurious match:', paste(spumatch), '\n')
    #print(pairdist[spumatch,])
  }
  cat('\n Total Number of mismatched samples =', length(nonmatch), '\n\n')
  
  
  
  ########### Attribute prediction and clinical samples swapping detection
  #traincli <- clinic[, c('sample', 'gender', 'msi')]
  #traincli$gender_prob <- apply(traincli, 1, function(x) if (x['gender'] == 'Female') 0 else 1)
  #traincli$msi_prob <- apply(traincli, 1, function(x) if (x['msi'] == 'MSI-High') 0 else 1)
  
  traincli <- clinic[, c('sample', clinical_attributes)]
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  for(i in 1:length(clinical_attributes)){
    #cli_attr <- clinical_attributes[i]
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    traincli[,cli_attr_prob_names_true[i]] <- as.factor(traincli[,cli_attr_prob_names_true[i]])
    traincli[,cli_attr_prob_names_true[i]] <- apply(traincli, 1, function(x) if (x[clinical_attributes[i]] == levels(traincli[,clinical_attributes[i]])[1]) 0 else 1)
  }
  
  #### First round of prediction using cross validation (after removing RNA/PRO mismatch samples)
  traincli <- predictCV(traincli, nonmatch, rna_sex, rna_atsm, pro_sex, pro_atsm, clinical_attributes)

  cat('\n')
  #### First round flagging potential Clinical Swapping
  #traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob) / 2
  #traincli$pred_msi    <- (traincli$rmsi_prob + traincli$pmsi_prob) / 2
  
  # prediction probability
  cli_attr_prob_names_pred_rna <- paste("r",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_pro <- paste("p",clinical_attributes,"_prob",sep = "")
  
  cli_attr_prob_names_pred_combine <- paste("pred_",clinical_attributes,sep = "")
  for(i in 1:length(clinical_attributes)){
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    traincli[,cli_attr_prob_names_pred_combine[i]] <- (traincli[,cli_attr_prob_names_pred_rna[i]] + traincli[,cli_attr_prob_names_pred_pro[i]]) / 2.0
  }
  
  
  
  clinic_swap <- flagClinicalSwap(traincli, nonmatch, clinical_attributes)
  cat('Clinical Label potential swapping:', paste0(clinic_swap), '\n\n')
  
  #### Second round of prediction (after removing RNA/PRO mismatched samples & clinical swapping cases)
  traincli <- predictLR(traincli, c(nonmatch, clinic_swap), rna_sex, rna_atsm, pro_sex, pro_atsm, clinical_attributes)
  
  write.table(traincli,file = paste(out_dir,"/clinical_attributes_pred.tsv",sep = ""),col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  
  ## output intermediate prediction file
  # for example: "sample","gender","GenderProvided","genderByRNA","genderByPRO","msi","MsiProvided","msiByRNA","msiByPRO"
  #tobewritten <- traincli[, c(1, 2, 4, 7, 11, 3, 5, 9, 13)]
  #colnames(tobewritten)[3:5] <- c('genderProvided', 'genderByRNA', 'genderByPRO')
  #colnames(tobewritten)[7:9] <- c('msiProvided', 'msiByRNA', 'msiByPRO')
  #intermediate_file <- paste(out_dir,"/intermediate.csv",sep = "")
  #write.table(tobewritten, intermediate_file, col.names=TRUE, row.names=FALSE, sep=',')
  

  final_tab <- data.frame(sample=sample_labels, Clinical=1:sample_n, RNAseq=1:sample_n, Proteomics=1:sample_n)
  
  cat('\n')
  ########### Determine Clinical Swapping (Second round to detect and correct the labels)
  #traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob)/2
  #traincli$pred_msi    <- (traincli$rmsi_prob + traincli$pmsi_prob) / 2
  for(i in 1:length(clinical_attributes)){
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    traincli[,cli_attr_prob_names_pred_combine[i]] <- (traincli[,cli_attr_prob_names_pred_rna[i]] + traincli[,cli_attr_prob_names_pred_pro[i]]) / 2.0
  }
  
  final_tab   <- correctClinicalSwapping(traincli, final_tab, nonmatch, clinical_attributes)
  clinic_swap <- which(final_tab$Clinical != 1:sample_n)
  
  cat('\n')
  ########## Determine omics swapping
  final_tab <- correctOmicsSwapping(traincli, final_tab, nonmatch, rnamatch, pairdist, clinical_attributes)
  rnaswap <- which(final_tab$RNAseq != 1:sample_n)
  proswap <- which(final_tab$Proteomics != 1:sample_n)
  swapped <- c(rnaswap, proswap)
  
  
  cat('\n')
  ########## Determine omics duplication and shifting
  dup_shift <- setdiff(nonmatch, swapped)
  final_tab <- correctOmicsShifting(traincli, final_tab, dup_shift, rankdist, rnamatch, promatch, swapped, pairdist, clinical_attributes)
  rnashift  <- which(final_tab$RNAseq != 1:sample_n)
  rnashift  <- setdiff(rnashift, rnaswap)
  proshift  <- which(final_tab$Proteomics != 1:sample_n)
  proshift  <- setdiff(proshift, proswap)
  
  
  ########## Output error statistics and corrected label
  errors <- c(length(clinic_swap), length(proswap), length(proshift), length(rnaswap), length(rnashift))
  names(errors) <- c('cli_swap', 'pro_swap', 'pro_shift', 'rna_swap', 'rna_shift')
  as.data.frame(errors) -> errors
  print(errors)

  cat('\n')
  final_tab_file <- paste0(out_dir, "/final.csv")
  write.table(final_tab, final_tab_file , col.names=TRUE, row.names=FALSE, sep=',')
  error_file <- paste0(out_dir, "/error.tsv")
  write.table(errors, error_file, row.names=T, col.names=T, sep='\t')
}




## Preprocessing: annotate gene with chromosome information
prpr_annotate <- function(geneSymbol, out_gene_file){
  cat('Annotating genes with chromosomes...\n')
  
  gene_info <- tryCatch({
    cat("Try www.ensembl.org ...\n")
    mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")     #ensemblRedirect = FALSE)
    gene_info_tmp <- getBM(attributes=c('hgnc_symbol', 'chromosome_name'),
                       filters='hgnc_symbol',
                       values=geneSymbol,
                       mart=mart)
    return(gene_info_tmp)
  },error=function(e){
    cat("Try http://uswest.ensembl.org/ ...\n")
    mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "http://uswest.ensembl.org/")     #ensemblRedirect = FALSE)
    gene_info_tmp <- getBM(attributes=c('hgnc_symbol', 'chromosome_name'),
                       filters='hgnc_symbol',
                       values=geneSymbol,
                       mart=mart)
    return(gene_info_tmp)
  },
  warning=function(cond){
    print(cond)
    message("Please see the file: kfoldCrossValidation.rda for related data!")
    save(e,ann_use,net_data,r,kfold,ranNum,file="kfoldCrossValidation.rda")
    return(NULL)
  })
  
  cat('Annotation is complete, save annotation into file ', out_gene_file, '\n\n')
  write.table(gene_info, out_gene_file, sep='\t', col.names=T, row.names=F)
}


## Preprocessing: proteomic
prpr_proteome <- function(proteome, sexgenes) {
  cat('Preprocessing proteomics... \n')
  
  ## partition proteome into sex genes and autosomal genes
  sex_pro  <- proteome[intersect(rownames(proteome), sexgenes), ]
  auto_pro <- proteome[setdiff(rownames(proteome), sexgenes), ]
  
  # Remove any sex gene which expression values is NA in all samples
  pmiss <- apply(sex_pro, 1, function(x) sum(is.na(x)))
  pmiss <- which(pmiss == ncol(sex_pro))
  if (length(pmiss) > 0) {
    cat('  ', length(pmiss), 'sex gene(s) has NA in ALL samples - removed... \n')
    sex_pro <- sex_pro[-pmiss, ]
  }
  # Replace missing value as 0 in sex genes
  sex_pro <- t(sex_pro)
  sex_pro[is.na(sex_pro)] <- 0
  
  # Remove any automosal gene with > 50% missing values
  pmiss <- apply(auto_pro, 1, function(x) sum(is.na(x)))
  pmiss <- which(pmiss > (ncol(auto_pro)/2))
  if (length(pmiss) > 0) {
    cat('  ', length(pmiss), 'of autosomal gene(s) has NA in > 50% samples - removed... \n')
    auto_pro <- auto_pro[-pmiss, ]
  }
  # Impute missing value in autosomal genes if >= 30% of genes has missing values; else remove
  missingrow <- nrow(auto_pro) - nrow(na.omit(auto_pro))
  missingPct <- missingrow / nrow(auto_pro)
  if (missingPct >= 0.30) {
    cat('  ', sprintf('%d(%.2f%%) of autosomal gene(s) has NA missing values - impute missing values... \n', missingrow, missingPct))
    auto_pro <- t(auto_pro)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    auto_imp <- missForest(auto_pro, parallelize="variables")
    stopCluster(cl)
    auto_pro <- auto_imp$ximp
  } else {
    cat('  ', sprintf('%d(%.2f%%) of autosomal gene(s) has NA missing values - removed... \n', missingrow, missingPct))
    auto_pro <- na.omit(auto_pro)
    auto_pro <- t(auto_pro)
  }
  
  proteome <- cbind(sex_pro, auto_pro)
  return(proteome)
}


## preprocessing: rnaseq
prpr_rnaseq <- function(rnaseq, sexgenes){
  cat('Preprocessing RNAseq... \n')
  
  ## partition into sex genes and autosomal genes
  sex_rna  <- rnaseq[intersect(rownames(rnaseq), sexgenes), ]
  auto_rna <- rnaseq[setdiff(rownames(rnaseq), sexgenes), ]
  
  # Remove any sex gene which expression values is NA in all samples
  pmiss <- apply(sex_rna, 1, function(x) sum(is.na(x)))
  pmiss <- which(pmiss == ncol(sex_rna))
  if (length(pmiss) > 0) {
    cat('  ', length(pmiss), 'sex gene(s) has NA in ALL samples - removed... \n')
    sex_rna <- sex_rna[-pmiss, ]
  }
  # Replace missing value as 0 in sex genes
  sex_rna <- t(sex_rna)
  sex_rna[is.na(sex_rna)] <- 0
  
  # Remove any automosal gene with > 50% missing values
  pmiss <- apply(auto_rna, 1, function(x) sum(is.na(x)))
  pmiss <- which(pmiss > (ncol(auto_rna)/2))
  if (length(pmiss) > 0) {
    cat('  ', length(pmiss), 'of autosomal gene(s) has NA in > 50% samples - removed... \n')
    auto_rna <- auto_rna[-pmiss, ]
  }
  # Impute missing value in autosomal genes if >= 30% of genes has missing values; else remove
  missingrow <- nrow(auto_rna) - nrow(na.omit(auto_rna))
  missingPct <- missingrow / nrow(auto_rna)
  if (missingPct >= 0.30) {
    cat('  ', sprintf('%d(%.2f%%) of autosomal gene(s) has NA missing values - impute missing values... \n', missingrow, missingPct))
    auto_rna <- t(auto_rna)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    auto_imp <- missForest(auto_rna, parallelize="variables")
    stopCluster(cl)
    auto_rna <- auto_imp$ximp
  } else {
    cat('  ', sprintf('%d(%.2f%%) of autosomal gene(s) has NA missing values - removed... \n', missingrow, missingPct))
    auto_rna <- na.omit(auto_rna)
    auto_rna <- t(auto_rna)
  }
  
  rnaseq <- cbind(sex_rna, auto_rna)
  return(rnaseq)
}


#### Preprocessing Module
preprocess <- function(rna_file, pro_file, gene_file, out_dir="./"){
  
  ## load protein and RNA expression data
  proteome <- read.delim(pro_file, stringsAsFactors = FALSE)
  rnaseq   <- read.delim(rna_file, stringsAsFactors = FALSE)
  
  ## annotate genes with chromosomes
  geneSymbol    <- union(rownames(proteome), rownames(rnaseq))
  out_gene_file <- paste0(out_dir, "/genes.tsv")
  gene_data <- read.delim(gene_file,stringsAsFactors = FALSE)
  write.table(gene_data,file = out_gene_file,row.names = FALSE,col.names = TRUE,sep = "\t",quote = FALSE)
  #prpr_annotate(geneSymbol, out_gene_file)
  sexgenes <- getSexGenes(out_gene_file)
  
  ## preprocessing rnaseq
  rnaseq   <- prpr_rnaseq(rnaseq, sexgenes)
  out_RNA  <- paste0(out_dir, "/rnaseq.tsv")
  write.table(rnaseq, out_RNA, row.names=TRUE, col.names=TRUE, sep='\t')
  
  ## preprocessing proteomic
  proteome <- prpr_proteome(proteome, sexgenes)
  out_PRO  <- paste0(out_dir, "/proteomic.tsv")
  write.table(proteome, out_PRO, row.names=TRUE, col.names=TRUE, sep='\t')
  
  return(c(out_RNA, out_PRO, out_gene_file))
}




## Load: gene chromosome annotation
getSexGenes <- function(gene_file){
  gene_info <- read.delim(gene_file, stringsAsFactors = FALSE)
  sexgenes <- which(gene_info$chromosome_name %in% c('X', 'Y'))
  sexgenes <- gene_info$hgnc_symbol[sexgenes]
  return(sexgenes)
}


## Load sample annotation data: test_cli.tsv
getClinical <- function(anno_file){
  clinic  <- read.delim(anno_file, stringsAsFactors = FALSE)
  return(clinic)
}


## Load RNAseq
getRNAseq <- function(rna_file, gene_file){
  rnaseq   <- read.delim(rna_file, stringsAsFactors = FALSE)
  sexgenes <- getSexGenes(gene_file)
  
  rna_atsm <- rnaseq[, setdiff(colnames(rnaseq), sexgenes)]
  rna_sex  <- rnaseq[, intersect(colnames(rnaseq), sexgenes)]
  
  rnaseq <- list()
  rnaseq$rna_atsm <- scale(rna_atsm)
  rnaseq$rna_sex  <- as.matrix(rna_sex)
  return(rnaseq)
}


## Load proteomic
getProteome <- function(pro_file, gene_file){
  proteome <- read.delim(pro_file, stringsAsFactors = FALSE)
  sexgenes <- getSexGenes(gene_file)
  
  pro_atsm <- proteome[, setdiff(colnames(proteome), sexgenes)]
  pro_sex  <- proteome[, intersect(colnames(proteome), sexgenes)]
  
  proteome <- list()
  proteome$pro_atsm <- scale(pro_atsm)
  proteome$pro_sex  <- as.matrix(pro_sex)
  return(proteome)
}




## Screen: compute correlation matrix between samples
computeCorr <- function(rna_atsm, pro_atsm, nonmatch = c()) {
  intergene <- intersect(colnames(rna_atsm), colnames(pro_atsm))
  if (length(nonmatch) == 0) {
    corgene <- extractCorGene(rna_atsm[, intergene], pro_atsm[, intergene])
  } else {
    corgene <- extractCorGene(rna_atsm[-nonmatch, intergene], pro_atsm[-nonmatch, intergene])
  }
  corsample <- cor(t(rna_atsm[, corgene]), t(pro_atsm[, corgene]))
  return(corsample)
}


## Screen: extract correlated genes between omics data
extractCorGene <- function(rnamatrix, promatrix){
  cormatrix <- cor(as.matrix(rnamatrix), as.matrix(promatrix))
  correlate <- diag(cormatrix)
  corgene <- names(which(correlate > 0.5))
  cat(sprintf('%d genes are highly correlated \n', length(corgene)))
  return(corgene)
}


## Screen: convert correlation matrix into probability matrix
computeRankdist <- function(corsample){
  rnadistpro <- corsample
  prodistrna <- corsample
  for (i in 1:nrow(corsample)){
    rnadistpro[i,] <- exp(scale(rnadistpro[i,])) / sum(exp(scale(rnadistpro[i,])))
    prodistrna[,i] <- exp(scale(prodistrna[,i])) / sum(exp(scale(prodistrna[,i])))
  }
  rankdist <- rnadistpro * prodistrna
  
  return(rankdist)
}


## Screen: Perform Stable Matching between omics samples and return the matching with scores
getMatching <- function(probmatrix){
  malerank <- rankbyRow(probmatrix)
  fmlerank <- rankbyCol(probmatrix)
  
  matcher <- stableMarriage(malerank, fmlerank)
  rnamatch <- as.numeric(sub('PRO_', '', matcher$malematch))
  promatch <- as.numeric(sub('RNA_', '', matcher$fmlematch))
  
  pairdist <- data.frame(rna=1:length(rnamatch), pro=rnamatch)
  pairdist$rnarank <- apply(pairdist, 1, function(x) malerank[x['rna'], x['pro']])
  pairdist$prorank <- apply(pairdist, 1, function(x) fmlerank[x['rna'], x['pro']])
  pairdist$distance <- pairdist$rnarank + pairdist$prorank
  pairdist$softmax <- apply(pairdist, 1, function(x) probmatrix[x['rna'], x['pro']])
  
  return(pairdist)
}


## StableMatching: get preferential rank by row
rankbyRow <- function(rankdist){
  malerank <- t(apply(rankdist, 1, function(x) rank(-x, ties.method='first')))
  colnames(malerank) <- colnames(rankdist)
  return(malerank)
}


## StableMatching: get preferential rank by column
rankbyCol <- function(rankdist) {
  fmlerank <- apply(rankdist, 2, function(x) rank(-x, ties.method='first'))
  rownames(fmlerank) <- rownames(rankdist)
  return(fmlerank)
}


## StableMatching: get matching pairs of samples with preferential ranks
stableMarriage <- function(malerank, fmlerank){
  matcher <- list()
  
  malematch <- rep(NA, nrow(malerank))
  names(malematch) <- rownames(malerank)
  fmlematch <- rep(NA, ncol(fmlerank))
  names(fmlematch) <- colnames(fmlerank)
  
  singlemales <- names(malematch)
  while (length(singlemales) != 0){
    for (ppsmale in singlemales){
      propose <- 1
      single <- TRUE
      while (single == TRUE){
        ppsfmle <- names(which(malerank[ppsmale,] == propose))
        engaged <- fmlematch[ppsfmle]
        if (is.na(engaged) || fmlerank[engaged, ppsfmle] > fmlerank[ppsmale, ppsfmle]) {
          malematch[ppsmale] <- ppsfmle
          fmlematch[ppsfmle] <- ppsmale
          singlemales <- setdiff(singlemales, ppsmale)
          if (!(is.na(engaged))) singlemales <- c(singlemales, engaged) 
          single <- FALSE
        } else {
          propose <- propose + 1
        }
      }
    }
  }
  
  matcher$malematch <- malematch
  matcher$fmlematch <- fmlematch
  return(matcher)
}




## Train: get weight of training instance, distributed equally by class
getClassWeight <- function(labels){
  weight <- rep( 1, length(labels) )
  labelclass <- unique(labels)
  weight[labels == labelclass[1]] <- length(labels)/ 2 / sum(labels == labelclass[1])
  weight[labels == labelclass[2]] <- length(labels)/ 2 / sum(labels == labelclass[2])
  return(weight)
}


## Train: weighted L1 regularized Logistic Regression, using CV to determine best lambda
trainGLM <- function(msiLabel, rnamatrix, alpha){
  weight <- getClassWeight(msiLabel)
  if (sum(weight) != length(msiLabel)){
    cat('Error: sum(classweight) does not equal to length(msiLabel)! \n')
  }   # should be nrow(trainset)
  
  # perform cross validation of elasticnet to determine optimum lambda
  cv.glm <- cv.glmnet(as.matrix(rnamatrix), msiLabel, family="binomial", weights=weight, alpha=alpha)
  (best_lambda <- cv.glm$lambda.1se)
  fit <- glmnet(as.matrix(rnamatrix), msiLabel, family="binomial", weights=weight, alpha=alpha, lambda=best_lambda)
  
  return(fit)
}


## Train: use trainGLM fn in k-fold
trainGLMcv <- function(msiLabel, inputmtx, alpha, k = 5){
  flds <- createFolds(msiLabel, k = k, list = TRUE, returnTrain = FALSE)
  predoutput <- data.frame(att=msiLabel, att_prob=0)
  
  for (f in 1:k) {
    testidx <- flds[[f]]
    numvar <- 0
    iter   <- 0
    cat('  Training Fold', f, '- Optimizing Model...')
    while (numvar < 4 && iter < 50){
      fit1 <- trainGLM(msiLabel[-testidx], inputmtx[-testidx, ], 0.3)
      numvar1 <- sum(coef(fit1) > 0)
      if (numvar1 >= numvar) {
        fit <- fit1
        numvar <- numvar1
      }
      iter <- iter + 1
    }
    cat('Done. \n')
    predoutput$att[testidx] <- predict(fit, as.matrix(inputmtx[testidx, ]), type='class')[, 1]
    predoutput$att_prob[testidx] <- predict(fit, as.matrix(inputmtx[testidx, ]), type='response')[, 1]
  }
  return(predoutput)
}

is_gender = function(x){
  gd <- grepl(pattern = "gender",x = x, ignore.case = TRUE)
  return(gd)
}


## Train: use trainGLMcv fn for each clinical attribute
# clinical_attributes = c("gender")
predictCV <- function(traincli, nonmatch, rna_sex, rna_atsm, pro_sex, pro_atsm, 
                      clinical_attributes){
  
  
  # prediction probability
  cli_attr_prob_names_pred_rna <- paste("r",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_pro <- paste("p",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_rna <- paste("r",clinical_attributes,sep = "")
  cli_attr_prob_names_pro <- paste("p",clinical_attributes,sep = "")
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  
  for(i in 1:length(clinical_attributes)){
    if(is_gender(clinical_attributes[i])){
      ## gender prediction
      cat("Predicting ", clinical_attributes[i] ," from RNA... \n")
      predoutput_rna <- trainGLMcv(traincli[,clinical_attributes[i]][-nonmatch], rna_sex[-nonmatch, ], 0.3)
      cat("Predicting ", clinical_attributes[i] ," from Protein... \n")
      predoutput_pro <- trainGLMcv(traincli[,clinical_attributes[i]][-nonmatch], pro_sex[-nonmatch, ], 0.3)
      
    }else{
      ## non-gender attribute prediction
      cat("Predicting ", clinical_attributes[i] ," from RNA... \n")
      predoutput_rna <- trainGLMcv(traincli[,clinical_attributes[i]][-nonmatch], rna_atsm[-nonmatch, ], 0.3)
      cat("Predicting ", clinical_attributes[i] ," from Protein... \n")
      predoutput_pro <- trainGLMcv(traincli[,clinical_attributes[i]][-nonmatch], pro_atsm[-nonmatch, ], 0.3)
    }
    
    ## RNA
    traincli[,cli_attr_prob_names_rna[i]] <- traincli[,clinical_attributes[i]]
    traincli[,cli_attr_prob_names_pred_rna[i]] <- traincli[,cli_attr_prob_names_true[i]]
    traincli[,cli_attr_prob_names_rna[i]][-nonmatch] <- predoutput_rna$att
    traincli[,cli_attr_prob_names_pred_rna[i]][-nonmatch] <- predoutput_rna$att_prob
    
    ## Protein
    traincli[,cli_attr_prob_names_pro[i]] <- traincli[,clinical_attributes[i]]
    traincli[,cli_attr_prob_names_pred_pro[i]] <- traincli[,cli_attr_prob_names_true[i]]
    traincli[,cli_attr_prob_names_pro[i]][-nonmatch] <- predoutput_pro$att
    traincli[,cli_attr_prob_names_pred_pro[i]][-nonmatch] <- predoutput_pro$att_prob
    
  }
  
  return(traincli)
}


find_gender_label = function(clinical_attributes = c("gender")){
  
  gender_label <- NA
  for(i in clinical_attributes){
    if(is_gender(i)){
      gender_label <- i
    }
  }
  return(gender_label)
}


## Correction: Detect potential clinical swapping
flagClinicalSwap <- function(traincli, nonmatch, clinical_attributes) {
  
  gender_label <- find_gender_label(clinical_attributes)
  if(!is.na(gender_label)){
    gender_prob <- paste(gender_label,"_prob",sep = "")
    pred_gender <- paste("pred_",gender_label,sep = "")
    #high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
    high_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.7)
    (high_suspect <- setdiff(high_suspect, nonmatch))
    #cli_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.45)
    cli_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.45)
    (cli_suspect <- setdiff(cli_suspect, nonmatch))
  }else{
    cli_suspect <- c()
  }
  
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_combine <- paste("pred_",clinical_attributes,sep = "")
  
  
  clinic_swap <- c()
  if (length(cli_suspect) <= 1){
    cat('No Clinical Swapping Cases Found! \n')
  } else  {
    #subset <- traincli[cli_suspect, c('gender_prob', 'msi_prob', 'pred_gender', 'pred_msi')]
    subset <- traincli[cli_suspect, c(cli_attr_prob_names_true,cli_attr_prob_names_pred_combine)]
    
    
    clidist <- matrix(nrow=length(cli_suspect), ncol=length(cli_suspect))
    rownames(clidist) <- cli_suspect
    colnames(clidist) <- cli_suspect
    for (male in as.character(cli_suspect)) {
      for (fmle in as.character(cli_suspect)) {
        #clidist[male, fmle] <- (1 - abs(subset[male, 'gender_prob'] - subset[fmle, 'pred_gender'])) + (1 - abs(subset[male, 'msi_prob'] - subset[fmle, 'pred_msi']))
        clidist[male, fmle] <- 0
        for(i in 1:length(clinical_attributes)){
          clidist[male, fmle] <- clidist[male, fmle] + (1 - abs(subset[male, cli_attr_prob_names_true[i]] - subset[fmle, cli_attr_prob_names_pred_combine[i]]))
        }
      }
    }
    
    climale <- rankbyRow(clidist)
    clifmle <- t(climale)
    
    clinic_match <- stableMarriage(climale, clifmle)
    for (eachmatch in names(clinic_match$malematch)){
      match1 <- as.numeric(eachmatch)
      match2 <- as.numeric(clinic_match$malematch[eachmatch])
      if (match1 %in% clinic_swap)    next
      if (match1 != match2) {
        clinic_swap <- c(clinic_swap, match1, match2)
      }
    }
  }
  clinic_swap <- union(clinic_swap, high_suspect)
  return(clinic_swap)
}


## Train: use trainGLM fn for each clinical attribute
# clinical_attributes = c("gender")
predictLR <- function(traincli, nonmatch, rna_sex, rna_atsm, pro_sex, pro_atsm,
                      clinical_attributes) {
  
  # prediction probability
  cli_attr_prob_names_pred_rna <- paste("r",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_pro <- paste("p",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_rna <- paste("r",clinical_attributes,sep = "")
  cli_attr_prob_names_pro <- paste("p",clinical_attributes,sep = "")
  
  
  for(i in 1:length(clinical_attributes)){
    if(is_gender(clinical_attributes[i])){
      numvar <- 0
      iter <- 0
      cat('Training RNA model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4){
        fit <- trainGLM(traincli[,clinical_attributes[i]][-nonmatch], rna_sex[-nonmatch, ], 0.3)
        numvar <- sum(coef(fit) > 0)
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      #traincli$rgender <- predict(fit, rna_sex, type='class')[, 1]
      traincli[,cli_attr_prob_names_rna[i]] <- predict(fit, rna_sex, type='class')[, 1]
      #traincli$rgender_prob <- predict(fit, rna_sex, type='response')[, 1]
      traincli[,cli_attr_prob_names_pred_rna[i]] <- predict(fit, rna_sex, type='response')[, 1]
      
      numvar <- 0
      iter   <- 0
      cat('Training Protein model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4 && iter < 50){
        fit1 <- trainGLM(traincli[,clinical_attributes[i]][-nonmatch], pro_sex[-nonmatch, ], 0.3)
        numvar1 <- sum(coef(fit1) > 0)
        if (numvar1 > numvar) {
          fit <- fit1
          numvar <- numvar1
        }
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      #traincli$pgender <- predict(fit, pro_sex, type='class')[, 1]
      traincli[,cli_attr_prob_names_pro[i]] <- predict(fit, pro_sex, type='class')[, 1]
      #traincli$pgender_prob <- predict(fit, pro_sex, type='response')[, 1]
      traincli[,cli_attr_prob_names_pred_pro[i]] <- predict(fit, pro_sex, type='response')[, 1]
      
      
    }else{
      ## non-gender attribute prediction
      numvar <- 0
      iter <- 0
      cat('Training RNA model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4){
        fit <- trainGLM(traincli[,clinical_attributes[i]][-nonmatch], rna_atsm[-nonmatch, ], 0.3)
        numvar <- sum(coef(fit) > 0)
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      traincli[,cli_attr_prob_names_rna[i]] <- predict(fit, rna_atsm, type='class')[, 1]
      traincli[,cli_attr_prob_names_pred_rna[i]] <- predict(fit, rna_atsm, type='response')[, 1]
      
      numvar <- 0
      iter   <- 0
      cat('Training Protein model for ',clinical_attributes[i],"\n",sep = "")
      while (numvar < 4 && iter < 50){
        fit1 <- trainGLM(traincli[,clinical_attributes[i]][-nonmatch], pro_atsm[-nonmatch, ], 0.3)
        numvar1 <- sum(coef(fit1) > 0)
        if (numvar1 > numvar) {
          fit <- fit1
          numvar <- numvar1
        }
        iter <- iter + 1
      }
      cat('Done - Regularized model with', numvar, 'variables. \n')
      traincli[,cli_attr_prob_names_pro[i]] <- predict(fit, pro_atsm, type='class')[, 1]
      traincli[,cli_attr_prob_names_pred_pro[i]] <- predict(fit, pro_atsm, type='response')[, 1]
      
    }
    
  }
  
  return(traincli)
}


## Correction: Detect clinical swapping and correct label
correctClinicalSwapping <- function(traincli, final_tab, nonmatch, clinical_attributes) {

  gender_label <- find_gender_label(clinical_attributes)
  if(!is.na(gender_label)){
    gender_prob <- paste(gender_label,"_prob",sep = "")
    pred_gender <- paste("pred_",gender_label,sep = "")
    #high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
    high_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.7)
    (high_suspect <- setdiff(high_suspect, nonmatch))
    if (length(high_suspect) > 0) {
      cat('Highly suspected Clinical swap case:', paste(high_suspect), '\n')
      final_tab$Clinical[high_suspect] <- -1
    }
    #cli_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.35)
    cli_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.35)
    cli_suspect <- setdiff(cli_suspect, nonmatch)
  }else{
    cli_suspect <- c()
  }
  
  
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_combine <- paste("pred_",clinical_attributes,sep = "")
  
  clinic_swap <- c()
  if (length(cli_suspect) <= 1){
    cat('No Clinical Swapping Cases Found! \n')
  } else  {
    #subset <- traincli[cli_suspect, c('gender_prob', 'msi_prob', 'pred_gender', 'pred_msi')]
    subset <- traincli[cli_suspect, c(cli_attr_prob_names_true,cli_attr_prob_names_pred_combine)]
    
    clidist <- matrix(nrow=length(cli_suspect), ncol=length(cli_suspect))
    rownames(clidist) <- cli_suspect
    colnames(clidist) <- cli_suspect
    for (male in as.character(cli_suspect)) {
      for (fmle in as.character(cli_suspect)) {
        #clidist[male, fmle] <- (1 - abs(subset[male, 'gender_prob'] - subset[fmle, 'pred_gender'])) + (1 - abs(subset[male, 'msi_prob'] - subset[fmle, 'pred_msi']))
        clidist[male, fmle] <- 0
        for(i in 1:length(clinical_attributes)){
          clidist[male, fmle] <- clidist[male, fmle] + (1 - abs(subset[male, cli_attr_prob_names_true[i]] - subset[fmle, cli_attr_prob_names_pred_combine[i]]))
        }
      }
    }
    
    climale <- rankbyRow(clidist)
    clifmle <- t(climale)
    
    clinic_match <- stableMarriage(climale, clifmle)
    for (eachmatch in names(clinic_match$malematch)){
      match1 <- as.numeric(eachmatch)
      match2 <- as.numeric(clinic_match$malematch[eachmatch])
      if (match1 %in% clinic_swap)    next
      if (match1 != match2) {
        final_tab[match1, 'Clinical'] <- match2
        final_tab[match2, 'Clinical'] <- match1
        clinic_swap <- c(clinic_swap, match1, match2)
        cat(sprintf('Clinical Label Swapping: %d <--> %d \n', match1, match2))
      }
    }
    print(traincli[clinic_swap, ])
  }

  return(final_tab)
}


## Correct: Detect omics swapping and correct label
correctOmicsSwapping <- function(traincli, final_tab, nonmatch, rnamatch, pairdist,
                                 clinical_attributes) {
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  # prediction probability
  cli_attr_prob_names_pred_rna <- paste("r",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_pro <- paste("p",clinical_attributes,"_prob",sep = "")
  
  swapped <- c()
  for (r in nonmatch){
    if (r %in% swapped)    next
    p <- rnamatch[r]
    if (rnamatch[p] == r && rnamatch[r] != r) {
      if (sum(pairdist[c(r, p), c('rnarank', 'prorank')]) > 8){
        cat(sprintf('Spurious Pair: %d <--> %d, distance = %f \n', r, p, sum(pairdist[c(r, p), c('rnarank', 'prorank')])))
        next
      }
      
      swapped <- c(swapped, r, p)
      #subset <- traincli[c(r,p), c('gender_prob', 'rgender_prob', 'pgender_prob', 'msi_prob', 'rmsi_prob', 'pmsi_prob')]
      #pro_swap <- abs(subset[1,1]-subset[1,2]) + abs(subset[1,1]-subset[2,3]) + abs(subset[2,1]-subset[2,2]) + abs(subset[2,1]-subset[1,3])
      #pro_swap <- pro_swap + abs(subset[1,4]-subset[1,5]) + abs(subset[1,4]-subset[2,6]) + abs(subset[2,4]-subset[2,5]) + abs(subset[2,4]-subset[1,6])
      #rna_swap <- abs(subset[1,1]-subset[2,2]) + abs(subset[1,1]-subset[1,3]) + abs(subset[2,1]-subset[1,2]) + abs(subset[2,1]-subset[2,3])
      #rna_swap <- rna_swap + abs(subset[1,4]-subset[2,5]) + abs(subset[1,4]-subset[1,6]) + abs(subset[2,4]-subset[1,5]) + abs(subset[2,4]-subset[2,6])
      
      pro_swap <- 0
      rna_swap <- 0
      for(i in 1:length(clinical_attributes)){
        subset <- traincli[c(r,p), c(cli_attr_prob_names_true[i], cli_attr_prob_names_pred_rna[i], cli_attr_prob_names_pred_pro[i])]
        pro_swap <- pro_swap + abs(subset[1,1]-subset[1,2]) + abs(subset[1,1]-subset[2,3]) + abs(subset[2,1]-subset[2,2]) + abs(subset[2,1]-subset[1,3])
        rna_swap <- rna_swap + abs(subset[1,1]-subset[2,2]) + abs(subset[1,1]-subset[1,3]) + abs(subset[2,1]-subset[1,2]) + abs(subset[2,1]-subset[2,3])
      }
      
      
      if (pro_swap < rna_swap){
        cat(sprintf('Proteome swap: %d <--> %d (RNA: %.3f vs PRO: %.3f) \n', r, p, rna_swap, pro_swap))
        final_tab[r, 'Proteomics'] <- p
        final_tab[p, 'Proteomics'] <- r
      } else {
        cat(sprintf('RNA swap: %d <--> %d (RNA: %.3f vs PRO: %.3f) \n', r, p, rna_swap, pro_swap))
        final_tab[r, 'RNAseq'] <- p
        final_tab[p, 'RNAseq'] <- r
      }
    }
  }
  return(final_tab)
}


## Correct: Detect omics duplication + shifting and correct label
correctOmicsShifting <- function(traincli, final_tab, dup_shift, rankdist, 
                                 rnamatch, promatch, swapped, pairdist,
                                 clinical_attributes) {
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  # prediction probability
  cli_attr_prob_names_pred_rna <- paste("r",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_pro <- paste("p",clinical_attributes,"_prob",sep = "")
  
  rnashift <- c()
  proshift <- c()
  if (length(dup_shift) == 0) {
    cat('No shifting & duplication suspect! Skip... \n')
  } else {
    cat(length(dup_shift), 'samples suspected for duplication and shifting:', paste(dup_shift), '\n')
    shiftdist <- pairdist[dup_shift,]
    shiftdist <- shiftdist[order(-shiftdist$distance), ]
    print(shiftdist)
    
    lose_starts <- shiftdist$rna[shiftdist$distance > 3]
    lose_ends   <- shiftdist$pro[shiftdist$distance > 3]
    
    
    ### chain identification
    chains <- list()
    i <- 1
    for (start in lose_starts){
      chain <- c(start)
      cnext <- start
      while (!(cnext %in% lose_ends)) {
        cnext <- promatch[cnext]
        chain <- c(chain, cnext)
      }
      
      chain <- c(chain, which.min(rankdist[, chain[length(chain)]]))
      chain <- c(which.min(rankdist[chain[1],]), chain)
      
      chains[[i]] <- chain
      i <- i + 1
    }
    
    if (sum(!(dup_shift %in% unlist(chains))) == 0){
      cat('All suspected samples are found in chain! \n\n')
    } else {
      cat('Warning! Not all suspected samples are found in chain. Circular shifting suspected! \n\n')
    }
    
    
    for (chain in chains){
      cat(paste(chain, collapse = ' --> '), '\n')
      if (chain[1] %in% swapped) {
        cat('Warning: chain head', chain[1], 'is swapped samples!\n')
      }
      if (chain[length(chain)] %in% swapped) {
        cat('Warning: chain tail', chain[length(chain)], 'is swapped samples!\n')
      }
      
      #subset <- traincli[chain, c('gender_prob', 'rgender_prob', 'pgender_prob', 'msi_prob', 'rmsi_prob', 'pmsi_prob')]
      #lenchain <- length(chain)
      #rna_shift <- sum(abs(subset[1:(lenchain-2), 1] - subset[2:(lenchain-1), 2])) + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 3])) + abs(subset[1,1] - subset[1,2])
      #rna_shift <- rna_shift + sum(abs(subset[1:(lenchain-2), 4] - subset[2:(lenchain-1), 5])) + sum(abs(subset[2:lenchain, 4] - subset[2:lenchain, 6])) + abs(subset[1,4] - subset[1,5])
      #pro_shift <- sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 2])) + sum(abs(subset[3:lenchain, 1] - subset[2:(lenchain-1), 3])) + abs(subset[lenchain,1] - subset[lenchain,3])
      #pro_shift <- pro_shift + sum(abs(subset[2:lenchain, 4] - subset[2:lenchain, 5])) + sum(abs(subset[3:lenchain, 4] - subset[2:(lenchain-1), 6])) + abs(subset[lenchain,4] - subset[lenchain,6])
      
      rna_shift <- 0
      pro_shift <- 0
      lenchain <- length(chain)
      for(i in 1:length(clinical_attributes)){
        subset <- traincli[chain, c(cli_attr_prob_names_true[i], cli_attr_prob_names_pred_rna[i], cli_attr_prob_names_pred_pro[i])]
        rna_shift <- rna_shift + sum(abs(subset[1:(lenchain-2), 1] - subset[2:(lenchain-1), 2])) + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 3])) + abs(subset[1,1] - subset[1,2])
        pro_shift <- pro_shift + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 2])) + sum(abs(subset[3:lenchain, 1] - subset[2:(lenchain-1), 3])) + abs(subset[lenchain,1] - subset[lenchain,3])
      }
      ## if rna_shift < pro_shift, means error rate for rna shifting is lower, means it is rna shifting
      
      distfront <- rankdist[chain[2], chain[1]]
      distback  <- rankdist[chain[lenchain], chain[lenchain-1]]
      ## if distback > distfront, means it is likely proteome duplication than RNAseq duplication
      
      mean_true <- c()
      for(i in 1:length(clinical_attributes)){
        subset <- traincli[chain, c(cli_attr_prob_names_true[i], cli_attr_prob_names_pred_rna[i], cli_attr_prob_names_pred_pro[i])]
        mean_true[i] <- mean(subset[, cli_attr_prob_names_true[i]][2:(lenchain-1)]) == 0 || mean(subset[, cli_attr_prob_names_true[i]][2:(lenchain-1)]) == 1
      }
      
      
      ## if all samples in a shifting chain the have same attribute (e.g. consistenly shifting ALL MALE & Low-MSI sample)
      #if ((mean(subset$gender_prob[2:(lenchain-1)]) == 0 || mean(subset$gender_prob[2:(lenchain-1)]) == 1) && (mean(subset$msi_prob[2:(lenchain-1)]) == 0 || mean(subset$msi_prob[2:(lenchain-1)]) == 1)) {
      if(all(mean_true)){
        if (distback > distfront) {
          cat('Same attr, so Proteome shift: ', chain[1], paste(chain[3:lenchain], collapse=' '), chain[length(chain)], sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: RNA_%d <-> PRO_%d = %.4f \t RNA_%d <-> PRO_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'Proteomics'] <- chain[3:length(chain)]
          proshift <- c(proshift, chain[2:(lenchain-1)])
        } else {
          cat('Same attr, so RNA shift: ', chain[1], paste(chain[1:(length(chain)-2)], collapse = ' '), chain[length(chain)], sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: RNA_%d <-> PRO_%d = %.4f \t RNA_%d <-> PRO_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'RNAseq'] <- chain[1:(length(chain)-2)]
          rnashift <- c(rnashift, chain[2:(lenchain-1)])
        }
        ## if all samples in a shifting chain the have different attribute (e.g. shifting of M/F or Low/high samples)
      } else{
        if (pro_shift < rna_shift) {          # supposedly distfront < distback
          cat('Proteome shift: ', chain[1], paste(chain[3:lenchain], collapse=' '), chain[length(chain)], sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: RNA_%d <-> PRO_%d = %.4f \t RNA_%d <-> PRO_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'Proteomics'] <- chain[3:length(chain)]
          proshift <- c(proshift, chain[2:(lenchain-1)])
          if (distfront > distback) {
            cat('Warning: distance indicates RNA_', chain[1], ' duplication but prediction results indicate Proteome shifting', sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift), '\n')
            heatmap(rankdist[chain, chain], Rowv=NA, Colv=NA)
          }
          
        } else if (rna_shift < pro_shift) {   # supposedly distfront > distback
          cat('RNA shift: ', chain[1], paste(chain[1:(length(chain)-2)], collapse = ' '), chain[length(chain)], sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: RNA_%d <-> PRO_%d = %.4f \t RNA_%d <-> PRO_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'RNAseq'] <- chain[1:(length(chain)-2)]
          rnashift <- c(rnashift, chain[2:(lenchain-1)])
          if (distfront < distback){
            cat('Warning: distance indicates PRO_', chain[lenchain], ' duplication but prediction results indicate RNAseq shifting', sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift), '\n')
            heatmap(rankdist[chain, chain], Rowv=NA, Colv=NA)
          }
        }
      }
      cat('\n')
    }
  }
  
  return(final_tab)
}


