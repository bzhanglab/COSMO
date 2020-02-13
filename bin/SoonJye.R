library(missForest)
library(biomaRt)
library(glmnet)
library(caret)
library(parallel)
library(doParallel)

run_2b = function(pro_file, rna_file, sample_file, out_dir="./"){
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ########### Data Retrieval
  cat('Correcting dataset:', pro_file, " and ", rna_file , '... \n')
  out_pro_files <- missing_value_impute(pro_file, rna_file, out_dir=out_dir)
  
  clinic  <- getClinical(sample_file)
  
  rnaseq <- getRNAseq(rna_file, out_pro_files[3])
  rna_atsm <- scale(rnaseq$rna_atsm)
  sample_n <- nrow(rna_atsm)
  cat('Total', sample_n, 'samples... \n')
  sample_labels <- rownames(rna_atsm)
  rownames(rna_atsm) <- paste0('RNA_', 1:sample_n)
  rna_sex  <- rnaseq$rna_sex
  rownames(rna_sex) <- paste0('RNA_', 1:sample_n)
  
  proteome <- getProteome(out_pro_files[1], out_pro_files[2])
  pro_atsm <- scale(proteome$pro_atsm)
  pro_sex  <- proteome$pro_sex
  
  
  ########### Screening RNA vs PRO mismatched samples with correlation
  intergene <- intersect(colnames(rna_atsm), colnames(pro_atsm))
  
  #### first round of getting correlated genes
  cat('Screening First Round: ')
  corgene <- extractCorGene(rna_atsm[, intergene], pro_atsm[, intergene])
  corsample <- cor(t(rna_atsm[, corgene]), t(pro_atsm[, corgene]))
  rankdist <- computeRankdist(corsample)
  pairdist <- getMatching(rankdist)
  
  nonmatch <- which(1:sample_n != pairdist$pro)
  cat('Screening First Round: Found', length(nonmatch), 'mismatched samples \n')
  
  #### second round of getting correlated genes (after removing nonmatch from first round)
  cat('Screening Second Round: \n')
  corgene <- extractCorGene(rna_atsm[-nonmatch, intergene], pro_atsm[-nonmatch, intergene])
  corsample <- cor(t(rna_atsm[, corgene]), t(pro_atsm[, corgene]))
  rankdist <- computeRankdist(corsample)
  pairdist <- getMatching(rankdist)
  
  rnamatch <- pairdist$pro
  promatch <- pairdist$rna[order(pairdist$pro)]
  nonmatch <- which(1:sample_n != rnamatch)
  cat('Screening Second Round: Found', length(nonmatch), 'mismatched samples \n')
  
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
    cat('Spurious match:', paste(spumatch), '\n')
    print(pairdist[spumatch,])
  }
  cat('\n Total Number of mismatched samples =', length(nonmatch), '\n\n')
  
  
  
  ########### Attribute prediction and clinical samples swapping detection
  traincli <- clinic[, c('sample', 'gender', 'msi')]
  traincli$gender_prob <- apply(traincli, 1, function(x) if (x['gender'] == 'Female') 0 else 1)
  traincli$msi_prob <- apply(traincli, 1, function(x) if (x['msi'] == 'MSI-High') 0 else 1)
  
  
  ########### First round of prediction (after removing RNA/PRO mismatch samples)
  predoutput <- trainGLMcv(traincli$gender[-nonmatch], rna_sex[-nonmatch, ], 0.3)
  traincli$rgender <- traincli$gender
  traincli$rgender_prob <- traincli$gender_prob
  traincli$rgender[-nonmatch] <- predoutput$att
  traincli$rgender_prob[-nonmatch] <- predoutput$att_prob
  
  
  predoutput <- trainGLMcv(traincli$msi[-nonmatch], rna_atsm[-nonmatch, ], 0.3)
  traincli$rmsi <- traincli$msi
  traincli$rmsi_prob <- traincli$msi_prob
  traincli$rmsi[-nonmatch] <- predoutput$att
  traincli$rmsi_prob[-nonmatch] <- predoutput$att_prob
  
  
  predoutput <- trainGLMcv(traincli$gender[-nonmatch], pro_sex[-nonmatch, ], 0.3)
  traincli$pgender <- traincli$gender
  traincli$pgender_prob <- traincli$gender_prob
  traincli$pgender[-nonmatch] <- predoutput$att
  traincli$pgender_prob[-nonmatch] <- predoutput$att_prob
  
  
  predoutput <- trainGLMcv(traincli$msi[-nonmatch], pro_atsm[-nonmatch, ], 0.3)
  traincli$pmsi <- traincli$msi
  traincli$pmsi_prob <- traincli$msi_prob
  traincli$pmsi[-nonmatch] <- predoutput$att
  traincli$pmsi_prob[-nonmatch] <- predoutput$att_prob
  
  
  cat('\n')
  ########### Determine Clinical Swapping (First round to identify potential swapping samples)
  traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob) / 2
  traincli$pred_msi <- (traincli$rmsi_prob + traincli$pmsi_prob) / 2
  
  high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
  (high_suspect <- setdiff(high_suspect, nonmatch))
  cli_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.45)
  (cli_suspect <- setdiff(cli_suspect, nonmatch))
  
  clinic_swap <- c()
  if (length(cli_suspect) <= 1){
    cat('No Clinical Swapping Cases Found! \n')
  } else  {
    cat('Clinical Label suspect swapping:', paste0(cli_suspect), '\n')
    subset <- traincli[cli_suspect, c('gender_prob', 'msi_prob', 'pred_gender', 'pred_msi')]
    
    clidist <- matrix(nrow=length(cli_suspect), ncol=length(cli_suspect))
    rownames(clidist) <- cli_suspect
    colnames(clidist) <- cli_suspect
    for (male in as.character(cli_suspect)) {
      for (fmle in as.character(cli_suspect)) {
        clidist[male, fmle] <- (1 - abs(subset[male, 'gender_prob'] - subset[fmle, 'pred_gender'])) + (1 - abs(subset[male, 'msi_prob'] - subset[fmle, 'pred_msi']))
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
  cat('Clinical Label potential swapping:', paste0(clinic_swap), '\n\n')
  
  
  ########### Second round of prediction (after removing RNA/PRO mismatched samples & clinical swapping cases)
  numvar <- 0
  iter <- 0
  cat('Training RNA_gender model... ')
  while (numvar < 4){
    fit <- trainGLM(traincli$gender[-c(nonmatch,clinic_swap)], rna_sex[-c(nonmatch,clinic_swap), ], 0.3)
    numvar <- sum(coef(fit) > 0)
    iter <- iter + 1
  }
  cat('Done - Regularized model with', numvar, 'variables. \n')
  traincli$rgender <- predict(fit, as.matrix(rna_sex), type='class')[, 1]
  traincli$rgender_prob <- predict(fit, as.matrix(rna_sex), type='response')[, 1]
  
  
  numvar <- 0
  iter <- 0
  cat('Training RNA_MSI model... ')
  while (numvar < 4){
    fit <- trainGLM(traincli$msi[-c(nonmatch,clinic_swap)], rna_atsm[-c(nonmatch,clinic_swap), ], 0.3)
    numvar <- sum(coef(fit) > 0)
    iter <- iter + 1
  }
  cat('Done - Regularized model with', numvar, 'variables. \n')
  traincli$rmsi <- predict(fit, as.matrix(rna_atsm), type='class')[, 1]
  traincli$rmsi_prob <- predict(fit, as.matrix(rna_atsm), type='response')[, 1]
  
  
  numvar <- 0
  iter   <- 0
  cat('Training PRO_gender model... ')
  while (numvar < 4 && iter < 50){
    fit1 <- trainGLM(traincli$gender[-c(nonmatch,clinic_swap)], pro_sex[-c(nonmatch,clinic_swap), ], 0.3)
    numvar1 <- sum(coef(fit1) > 0)
    if (numvar1 > numvar) {
      fit <- fit1
      numvar <- numvar1
    }
    iter <- iter + 1
  }
  cat('Done - Regularized model with', numvar, 'variables. \n')
  traincli$pgender <- predict(fit, as.matrix(pro_sex), type='class')[, 1]
  traincli$pgender_prob <- predict(fit, as.matrix(pro_sex), type='response')[, 1]
  
  
  numvar <- 0
  iter   <- 0
  cat('Training PRO_MSI model... ')
  while (numvar < 4 && iter < 50){
    fit1 <- trainGLM(traincli$msi[-c(nonmatch,clinic_swap)], pro_atsm[-c(nonmatch,clinic_swap), ], 0.3)
    numvar1 <- sum(coef(fit1) > 0)
    if (numvar1 > numvar) {
      fit <- fit1
      numvar <- numvar1
    }
    iter <- iter + 1
  }
  cat('Done - Regularized model with', numvar, 'variables. \n')
  traincli$pmsi <- predict(fit, as.matrix(pro_atsm), type='class')[, 1]
  traincli$pmsi_prob <- predict(fit, as.matrix(pro_atsm), type='response')[, 1]
  
  
  # output intermediate prediction file
  tobewritten <- traincli[, c(1, 2, 4, 7, 11, 3, 5, 9, 13)]
  colnames(tobewritten)[3:5] <- c('GenderProvided', 'GenderByRNA', 'GenderByPRO')
  colnames(tobewritten)[7:9] <- c('MsiProvided', 'MsiByRNA', 'MsiByPRO')
  intermediate_file <- paste(out_dir,"/intermediate.csv",sep = "")
  write.table(tobewritten, intermediate_file, col.names=TRUE, row.names=FALSE, sep=',')
  
  
  ########### Label Correction
  final_tab <- data.frame(sample=sample_labels, Clinical=1:sample_n, RNAseq=1:sample_n, Proteomics=1:sample_n)
  
  
  cat('\n')
  ########### Determine Clinical Swapping (Second round to detect again and correct the labels)
  traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob)/2
  traincli$pred_msi <- (traincli$rmsi_prob + traincli$pmsi_prob) / 2
  
  high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
  (high_suspect <- setdiff(high_suspect, nonmatch))
  if (length(high_suspect) > 0) {
    cat('Highly suspected Clinical swap case:', paste(high_suspect), '\n')
    final_tab$Clinical[high_suspect] <- -1
  }
  cli_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.35)
  cli_suspect <- setdiff(cli_suspect, nonmatch)
  
  
  clinic_swap <- c()
  if (length(cli_suspect) <= 1){
    cat('No Clinical Swapping Cases Found! \n')
  } else  {
    cat('Suspected Clinical swap case:', paste(cli_suspect), '\n')
    subset <- traincli[cli_suspect, c('gender_prob', 'msi_prob', 'pred_gender', 'pred_msi')]
    
    clidist <- matrix(nrow=length(cli_suspect), ncol=length(cli_suspect))
    rownames(clidist) <- cli_suspect
    colnames(clidist) <- cli_suspect
    for (male in as.character(cli_suspect)) {
      for (fmle in as.character(cli_suspect)) {
        clidist[male, fmle] <- (1 - abs(subset[male, 'gender_prob'] - subset[fmle, 'pred_gender'])) + (1 - abs(subset[male, 'msi_prob'] - subset[fmle, 'pred_msi']))
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
  
  clinic_swap <- union(clinic_swap, high_suspect)
  
  
  cat('\n')
  ########## Determine RNA or PRO swapping
  proswap <- c()
  rnaswap <- c()
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
      subset <- traincli[c(r,p), c('gender_prob', 'rgender_prob', 'pgender_prob', 'msi_prob', 'rmsi_prob', 'pmsi_prob')]
      pro_swap <- abs(subset[1,1]-subset[1,2]) + abs(subset[1,1]-subset[2,3]) + abs(subset[2,1]-subset[2,2]) + abs(subset[2,1]-subset[1,3])
      pro_swap <- pro_swap + abs(subset[1,4]-subset[1,5]) + abs(subset[1,4]-subset[2,6]) + abs(subset[2,4]-subset[2,5]) + abs(subset[2,4]-subset[1,6])
      rna_swap <- abs(subset[1,1]-subset[2,2]) + abs(subset[1,1]-subset[1,3]) + abs(subset[2,1]-subset[1,2]) + abs(subset[2,1]-subset[2,3])
      rna_swap <- rna_swap + abs(subset[1,4]-subset[2,5]) + abs(subset[1,4]-subset[1,6]) + abs(subset[2,4]-subset[1,5]) + abs(subset[2,4]-subset[2,6])
      
      if (pro_swap < rna_swap){
        cat(sprintf('Proteome swap: %d <--> %d (RNA: %.3f vs PRO: %.3f) \n', r, p, rna_swap, pro_swap))
        final_tab[r, 'Proteomics'] <- p
        final_tab[p, 'Proteomics'] <- r
        proswap <- c(proswap, r, p)
      } else {
        cat(sprintf('RNA swap: %d <--> %d (RNA: %.3f vs PRO: %.3f) \n', r, p, rna_swap, pro_swap))
        final_tab[r, 'RNAseq'] <- p
        final_tab[p, 'RNAseq'] <- r
        rnaswap <- c(rnaswap, r, p)
      }
    }
  }
  
  
  
  cat('\n')
  ########## Determine duplication and shifting
  dup_shift <- setdiff(nonmatch, swapped)
  rnashift <- c()
  proshift <- c()
  if (length(dup_shift) == 0) {
    cat('No shifting & duplication suspect! Skip... \n')
  } else {
    cat(length(dup_shift), 'samples suspected for duplication and shifting =', paste(dup_shift), '\n')
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
      
      subset <- traincli[chain, c('gender_prob', 'rgender_prob', 'pgender_prob', 'msi_prob', 'rmsi_prob', 'pmsi_prob')]
      lenchain <- length(chain)
      rna_shift <- sum(abs(subset[1:(lenchain-2), 1] - subset[2:(lenchain-1), 2])) + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 3])) + abs(subset[1,1] - subset[1,2])
      rna_shift <- rna_shift + sum(abs(subset[1:(lenchain-2), 4] - subset[2:(lenchain-1), 5])) + sum(abs(subset[2:lenchain, 4] - subset[2:lenchain, 6])) + abs(subset[1,4] - subset[1,5])
      pro_shift <- sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 2])) + sum(abs(subset[3:lenchain, 1] - subset[2:(lenchain-1), 3])) + abs(subset[lenchain,1] - subset[lenchain,3])
      pro_shift <- pro_shift + sum(abs(subset[2:lenchain, 4] - subset[2:lenchain, 5])) + sum(abs(subset[3:lenchain, 4] - subset[2:(lenchain-1), 6])) + abs(subset[lenchain,4] - subset[lenchain,6])
      ## if rna_shift < pro_shift, means error rate for rna shifting is lower, means it is rna shifting
      
      distfront <- rankdist[chain[2], chain[1]]
      distback  <- rankdist[chain[lenchain], chain[lenchain-1]]
      ## if distback > distfront, means it is likely proteome duplication than RNAseq duplication
      
      
      ## if all samples in a shifting chain the have same attribute (e.g. consistenly shifting ALL MALE & Low-MSI sample)
      if ((mean(subset$gender_prob[2:(lenchain-1)]) == 0 || mean(subset$gender_prob[2:(lenchain-1)]) == 1) && (mean(subset$msi_prob[2:(lenchain-1)]) == 0 || mean(subset$msi_prob[2:(lenchain-1)]) == 1)) {
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


## Annotate gene with chromosome information
prpr_annotate <- function(geneSymbol, out_gene_file){
  cat('Annotating genes with chromosomes...\n')
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")     #ensemblRedirect = FALSE)
  gene_info <- getBM(attributes=c('hgnc_symbol', 'chromosome_name'),
                     filters='hgnc_symbol',
                     values=geneSyms,
                     mart=mart)
  
  cat('Annotation is complete, save annotation into file ', out_gene_file, '\n')
  write.table(gene_info, out_gene_file, sep='\t', col.names=T, row.names=F)
}

#### Preprocessing module
missing_value_impute <- function(pro_file, rna_file, out_dir="./", gene_file=NULL){
  
  ## load protein and RNA expression data    
  proteome <- read.delim(pro_file, stringsAsFactors = FALSE)
  rnaseq   <- read.delim(rna_file, stringsAsFactors = FALSE)
  
  ## Annotate genes with chromosomes
  geneSyms <- union(rownames(proteome), rownames(rnaseq))
  out_gene_file <- paste0(out_dir, "/genes.tsv")
  prpr_annotate(geneSyms, out_gene_file)

  cat('Handling missing values of proteomics from', pro_file, "\n")
  
  ## extract sex genes to partition the datasets
  sexgenes <- getSexGenes(out_gene_file)
  ## partition proteome into sex genes and autosomal genes
  # replace missing value as 0 in sex genes
  sexprot <- proteome[intersect(rownames(proteome), sexgenes), ]
  sexprot <- t(sexprot)
  sexprot[is.na(sexprot)] <- 0
  out_f1 <- paste0(out_dir, "/sexprot.tsv")
  write.table(sexprot, out_f1, row.names=TRUE, col.names=TRUE, sep='\t')
  
  # remove genes with > 50% missing values
  autopro <- proteome[setdiff(rownames(proteome), sexgenes), ]
  pmiss <- apply(autopro, 1, function(x) sum(is.na(x)))
  pmiss <- which(pmiss < (nrow(autopro)/2))
  autopro <- autopro[pmiss, ]
  # impute remaining missing value in autosomal genes matrix
  autopro <- t(autopro)
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  autoimp <- missForest(autopro, parallelize="variables")
  stopCluster(cl)
  out_f2 <- paste0(out_dir, "/autoprot.tsv")
  write.table(autoimp$ximp, out_f2, row.names=TRUE, col.names=TRUE, sep='\t')
  
  return(c(out_f1, out_f2, out_gene_file))
}


getSexGenes <- function(gene_file){
  gene_info <- read.table(gene_file, header=T)
  sexgenes <- which(gene_info$chromosome_name %in% c('X', 'Y'))
  sexgenes <- gene_info$hgnc_symbol[sexgenes]
  return(sexgenes)
}


## read sample annotation data: test_cli.tsv
getClinical <- function(clifile){
  clinic  <- read.table(clifile, header=TRUE)
  return(clinic)
}


getRNAseq <- function(rna_file, gene_file){
  rnadata  <- read.delim(rna_file, stringsAsFactors = FALSE)
  sexgenes <- getSexGenes(gene_file)
  
  rna_atsm <- rnadata[setdiff(rownames(rnadata), sexgenes), ]
  rna_atsm <- na.exclude(rna_atsm)
  rna_atsm <- t(rna_atsm)
  
  rna_sex  <- rnadata[intersect(rownames(rnadata), sexgenes), ]
  rna_sex[is.na(rna_sex)] <- 0
  rna_sex  <- rna_sex[rowSums(rna_sex) != 0, ]
  rna_sex  <- t(rna_sex)
  
  rnaseq <- list()
  rnaseq$rna_atsm <- rna_atsm
  rnaseq$rna_sex  <- rna_sex
  return(rnaseq)
}


getProteome <- function(s_file, ns_file){
  profile <- ns_file
  pro_atsm <- read.table(profile, header=T)
  rownames(pro_atsm) <- paste0('PRO_', 1:nrow(pro_atsm))
  
  profile <- s_file
  pro_sex  <- read.table(profile, header=T)
  rownames(pro_sex) <- paste0('PRO_', 1:nrow(pro_sex))
  pro_sex  <- pro_sex[, colSums(pro_sex) != 0]
  
  proteome <- list()
  proteome$pro_atsm <- as.matrix(pro_atsm)
  proteome$pro_sex  <- as.matrix(pro_sex)
  return(proteome)
}


getClassWeight <- function(labels){
  weight <- rep( 1, length(labels) )
  labelclass <- unique(labels)
  weight[labels == labelclass[1]] <- length(labels)/ 2 / sum(labels == labelclass[1])
  weight[labels == labelclass[2]] <- length(labels)/ 2 / sum(labels == labelclass[2])
  return(weight)
}



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


## Perform stable matching between RNA and PRO samples and return matching score
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
  pairdist$softmax <- apply(pairdist, 1, function(x) rankdist[x['rna'], x['pro']])
  
  return(pairdist)
}

rankbyRow <- function(rankdist){
  malerank <- t(apply(rankdist, 1, function(x) rank(-x, ties.method='first')))
  colnames(malerank) <- colnames(rankdist)
  return(malerank)
}


rankbyCol <- function(rankdist) {
  fmlerank <- apply(rankdist, 2, function(x) rank(-x, ties.method='first'))
  rownames(fmlerank) <- rownames(rankdist)
  return(fmlerank)
}

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



extractCorGene <- function(rnamatrix, promatrix){
  cormatrix <- cor(as.matrix(rnamatrix), as.matrix(promatrix) )
  correlate <- diag(cormatrix)
  corgene <- names(which(correlate > 0.5))
  cat(sprintf('%d genes are highly correlated \n', length(corgene)))
  return(corgene)
}

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



trainGLMcv <- function(msiLabel, inputmtx, alpha, k = 5){
  flds <- createFolds(msiLabel, k = k, list = TRUE, returnTrain = FALSE)
  predoutput <- data.frame(att=msiLabel, att_prob=0)
  
  for (f in 1:k) {
    testidx <- flds[[f]]
    numvar <- 0
    iter   <- 0
    cat('Training Fold', f, '- Optimizing Model...')
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

cmdargs <- commandArgs(TRUE)
pro_file <- cmdargs[1]
rna_file <- cmdargs[2]
sample_file <- cmdargs[3]
out_dir <- cmdargs[4]
#prefix <- cmdargs[5]


run_2b(pro_file, rna_file, sample_file, out_dir=out_dir)


