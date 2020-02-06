library(missForest)
library(biomaRt)
library(glmnet)
library(caret)
library(parallel)
library("doParallel")


run_2c = function(pro_file, rna_file, sample_file, out_dir="./"){
    
    dir.create(out_dir,recursive = TRUE,showWarnings = FALSE)
    
    ## errors variable store the values of each type of error in each dataset
    ## errors <- matrix(0, nrow=50, ncol=6)
    ## colnames(errors) <- c('dataset', 'cli_swap', 'pro_swap', 'pro_shift', 'rna_swap', 'rna_shift')
    ## errors[, 1] <- 1:50
    
    ########### Data Preprocessing and Retrieval
    cat('Correcting dataset:', pro_file, " and ", rna_file , '... \n')
    out_pro_files <- missing_value_impute(pro_file,rna_file,out_dir=out_dir)
    
    clinic  <- getClinical(sample_file)
    
    colnames(clinic)[3] <- 'gender'
    
    rnaseq <- getRNAseq(rna_file,out_pro_files[3])
    
    rna_atsm <- scale(rnaseq$rna_atsm)
    rna_sex  <- rnaseq$rna_sex
    
    proteome <- getProteome(out_pro_files[1],out_pro_files[2])
    pro_atsm <- scale(proteome$pro_atsm)
    pro_sex  <- proteome$pro_sex
    
    
    sample_n <- nrow(rna_atsm)
    
    
    ########### Screening RNA vs PRO mismatch with correlation
    intergene <- intersect(colnames(rna_atsm), colnames(pro_atsm))
    
    
    #### first round of getting correlated genes
    corgene <- extractCorGene(rna_atsm[, intergene], pro_atsm[, intergene])
    corsample <- cor(t(rna_atsm[, corgene]), t(pro_atsm[, corgene]))
    probmatrix <- computeRankdist(corsample)
    pairdist <- getMatching(probmatrix)
    
    nonmatch <- which(1:80 != pairdist$pro)
    length(nonmatch)
    
    
    #### second round of getting correlated genes (after removing mismatch from first round)
    corgene <- extractCorGene(rna_atsm[-nonmatch, intergene], pro_atsm[-nonmatch, intergene])
    corsample <- cor(t(rna_atsm[, corgene]), t(pro_atsm[, corgene]))
    probmatrix <- computeRankdist(corsample)
    pairdist <- getMatching(probmatrix)
    
    nonmatch <- which(1:80 != pairdist$pro)
    length(nonmatch)
    
    
    malerank <- rankbyRow(probmatrix)
    fmlerank <- rankbyCol(probmatrix)
    #### Determine spurious match due to being left out (threshold is log(n))
    if (sum(1:80 == pairdist$pro & pairdist$distance >= log(sample_n)) > 0){
        spumatch <- pairdist$rna[1:80 == pairdist$pro & pairdist$distance >= log(80)]
        nonmatch <- c(nonmatch, spumatch)
        cat('Spurious match:', paste(spumatch), '\n')
        print(pairdist[spumatch,])
    }
    cat('\n Number of Non-matching =', length(nonmatch), '\n\n')
    
    
    #### output correlation file
    #write.table(corsample, sprintf('outdir/for_2c/Testing_%d/2c_correlation_%d.csv', d, d) , col.names=TRUE, row.names=TRUE, sep=',')
    
    
    
    
    ########### Attribute prediction and clinical samples swapping detection
    traincli <- clinic[, c('sample', 'gender')]
    traincli$gender_prob <- apply(traincli, 1, function(x) if (x['gender'] == 'Female') 0 else 1)
    
    
    #### First round of prediction (after removing RNA/PRO mismatch samples)
    predoutput <- trainGLMcv(traincli$gender[-nonmatch], rna_sex[-nonmatch, ], 0.3)
    traincli$rgender <- traincli$gender
    traincli$rgender_prob <- traincli$gender_prob
    traincli$rgender[-nonmatch] <- predoutput$gender
    traincli$rgender_prob[-nonmatch] <- predoutput$gender_prob
    
    
    predoutput <- trainGLMcv(traincli$gender[-nonmatch], pro_sex[-nonmatch, ], 0.3)
    traincli$pgender <- traincli$gender
    traincli$pgender_prob <- traincli$gender_prob
    traincli$pgender[-nonmatch] <- predoutput$gender
    traincli$pgender_prob[-nonmatch] <- predoutput$gender_prob
    
    
    
    cat('\n')
    #### Determine Clinical Swapping (First round to identify potential swapping samples)
    traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob) / 2
    
    high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
    (high_suspect <- setdiff(high_suspect, nonmatch))
    cli_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.45)
    cli_suspect <- setdiff(cli_suspect, nonmatch)
    
    clinic_swap <- determineClinicalSwapping(traincli, cli_suspect)
    clinic_swap <- union(clinic_swap, high_suspect)
    cat('1st iter: Potential Clinical Label swapping:', paste0(clinic_swap), '\n\n')
    
    
    
    #### Second round of prediction (after removing RNA/PRO mismatched samples & clinical swapping cases)
    numvar <- 0
    iter <- 0
    while (numvar < 4){
        fit <- trainGLM(traincli$gender[-c(nonmatch,clinic_swap)], rna_sex[-c(nonmatch,clinic_swap), ], 0.3)
        numvar <- sum(coef(fit) > 0)
        iter <- iter + 1
        cat('Iter', iter, '- RNAseq Gender predictor model has', numvar, 'variable. \n')
    }
    traincli$rgender <- predict(fit, as.matrix(rna_sex), type='class')[, 1]
    traincli$rgender_prob <- predict(fit, as.matrix(rna_sex), type='response')[, 1]
    
    
    numvar <- 0
    iter   <- 0
    while (numvar < 4 && iter < 50){
        fit1 <- trainGLM(traincli$gender[-c(nonmatch,clinic_swap)], pro_sex[-c(nonmatch,clinic_swap), ], 0.3)
        numvar1 <- sum(coef(fit1) > 0)
        if (numvar1 > numvar) {
            fit <- fit1
            numvar <- numvar1
        }
        cat('Iter', iter, ' - Proteomics Gender predictor model has', numvar, 'variable. \n')
        iter <- iter + 1
    }
    traincli$pgender <- predict(fit, as.matrix(pro_sex), type='class')[, 1]
    traincli$pgender_prob <- predict(fit, as.matrix(pro_sex), type='response')[, 1]
    
    
    #### output prediction intermediate file
    tobewritten <- traincli[, c(1, 2, 3, 5, 7)]
    colnames(tobewritten)[3:5] <- c('GenderProvided', 'GenderByRNA', 'GenderByPRO')
    #write.table(tobewritten, sprintf('outdir/for_2c/Testing_%d/2c_intermediate_%d.csv', d, d) , col.names=TRUE, row.names=FALSE, sep=',')
    
    
    ########### Label Correction
    final_tab <- data.frame(sample=paste0('Testing_', 1:sample_n), Clinical=1:sample_n, RNAseq=1:sample_n, Proteomics=1:sample_n)
    
    cat('\n')
    #### Determine Clinical Swapping (Second round to detect again and correct the labels)
    traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob)/2
    
    high_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.7)
    (high_suspect <- setdiff(high_suspect, nonmatch))
    if (length(high_suspect) > 0) {
        cat('Highly suspected Clinical swap case:', paste(high_suspect), '\n')
        final_tab$Clinical[high_suspect] <- -1
    }
    cli_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.35)
    cli_suspect <- setdiff(cli_suspect, nonmatch)
    
    
    clinic_swap <- determineClinicalSwapping(traincli, cli_suspect)
    for (i in seq(1, length(clinic_swap), 2)){
        r <- clinic_swap[i]
        p <- clinic_swap[i+1]
        final_tab$Clinical[r] <- p
        final_tab$Clinical[p] <- r
    }
    clinic_swap <- union(clinic_swap, high_suspect)
    cat('2nd iter: Final Clinical Label swapping:', paste0(clinic_swap), '\n\n')
    
    
    cat('\n')
    #### Determine RNA or PRO swapping
    proswap <- c()
    rnaswap <- c()
    swapped <- c()
    for (r in nonmatch){
        if (r %in% swapped)    next
        p <- pairdist$pro[r]
        if (pairdist$pro[p] == r && pairdist$pro[r] != r) {
            if (malerank[r,p] + fmlerank[r,p] + malerank[p,r] + fmlerank[p,r] > 8){
                cat(sprintf('Spurious Pair: %d <--> %d, distance = %f \n', r, p, malerank[r,p] + fmlerank[r,p] + malerank[p,r] + fmlerank[p,r]))
                next
            }
            
            swapped <- c(swapped, r, p)
            subset <- traincli[c(r,p), c('gender_prob', 'rgender_prob', 'pgender_prob')]
            
            # calculate error rate of swapping PRO and error rate of swapping RNA
            pro_swap <- abs(subset[1,1]-subset[1,2]) + abs(subset[1,1]-subset[2,3]) + abs(subset[2,1]-subset[2,2]) + abs(subset[2,1]-subset[1,3])
            rna_swap <- abs(subset[1,1]-subset[2,2]) + abs(subset[1,1]-subset[1,3]) + abs(subset[2,1]-subset[1,2]) + abs(subset[2,1]-subset[2,3])
            # if pro_swap < rna_swap, means error rate is lower swapping PRO than swapping RNA, hence it is likely PRO swapping
            
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
    #### Determine duplication and shifting
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
        
        
        #### chain identification
        chains <- list()
        i <- 1
        for (start in lose_starts){
            chain <- c(start)
            cnext <- start
            while (!(cnext %in% lose_ends)) {
                cnext <- pairdist$rna[pairdist$pro == cnext]
                chain <- c(chain, cnext)
            }
            
            chain <- c(chain, which.min(fmlerank[, chain[length(chain)]]))
            chain <- c(which.min(malerank[chain[1],]), chain)
            
            chains[[i]] <- chain
            i <- i + 1
        }
        
        if (sum(!(dup_shift %in% unlist(chains))) == 0){
            cat('All suspected samples are found in chain! \n\n')
        } else {
            cat('Warning! Not all suspected samples are found in chain. Circular shifting suspected! \n\n')
        }
        
        #### Determine if it is PRO or RNA chain shifting
        for (chain in chains){
            cat(paste(chain, collapse = ' --> '), '\n')
            if (chain[1] %in% swapped) {
                cat('Warning: chain head', chain[1], 'is swapped samples!\n')
            }
            if (chain[length(chain)] %in% swapped) {
                cat('Warning: chain tail', chain[length(chain)], 'is swapped samples!\n')
            }
            
            subset <- traincli[chain, c('gender_prob', 'rgender_prob', 'pgender_prob')]
            lenchain <- length(chain)
            
            rna_shift <- sum(abs(subset[1:(lenchain-2), 1] - subset[2:(lenchain-1), 2])) + sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 3])) + abs(subset[1,1] - subset[1,2])
            pro_shift <- sum(abs(subset[2:lenchain, 1] - subset[2:lenchain, 2])) + sum(abs(subset[3:lenchain, 1] - subset[2:(lenchain-1), 3])) + abs(subset[lenchain,1] - subset[lenchain,3])
            ## if rna_shift < pro_shift, means error rate for RNA shifting is lower, means it is RNA shifting
            
            distfront <- rankdist[chain[2], chain[1]]
            distback  <- rankdist[chain[lenchain], chain[lenchain-1]]
            ## if distback > distfront, means it is likely proteome duplication than RNAseq duplication
            
            
            ## if all samples in a shifting chain the have same attribute (e.g. consistenly shifting ALL MALE sample)
            if (mean(subset$gender_prob[2:(lenchain-1)]) == 0 || mean(subset$gender_prob[2:(lenchain-1)]) == 1) {
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
                ## if all samples in a shifting chain have different attribute (e.g. shifting of MALE & FEMALE sample)
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
    colnames(errors) <- c('dataset', 'cli_swap', 'pro_swap', 'pro_shift', 'rna_swap', 'rna_shift')
    as.data.frame(errors) -> errors
    print(errors)
    cat('\n')
    #write.table(final_tab, sprintf('outdir/for_2c/Testing_%d/final2cv4_%d.csv', d, d) , col.names=TRUE, row.names=FALSE, sep=',')
    #write.table(errors, 'outdir/for_2c/error2cv4.tsv', row.names=F, col.names=T, sep='\t')
    
    #write.table(final_tab, sprintf('data/for_2b/Testing_%d/final2bv2_%d.csv', d, d) , col.names=TRUE, row.names=FALSE, sep=',')
    final_tab_file <- paste(out_dir,"/final.csv",sep = "")
    write.table(final_tab, final_tab_file , col.names=TRUE, row.names=FALSE, sep=',')
    error_file <- paste(out_dir,"/error.tsv",sep = "")
    #write.table(errors, 'data/for_2b/error2bv2.tsv', row.names=F, col.names=T, sep='\t')
    write.table(errors, error_file, row.names=F, col.names=T, sep='\t')
    
}

## Impute missing value of proteomics
missing_value_impute=function(pro_file,rna_file,out_dir="./",gene_file=NULL){
    
    #proteome <- read.table(sprintf('../mislabeling/data/for_2c/Testing_%d/test_pro.tsv', d), header=T)
    # 11355 x 80
    #rnaseq   <- read.table(sprintf('../mislabeling/data/for_2c/Testing_%d/test_rna.tsv', d), header=T)
    # 19275 x 80
    
    proteome <- read.delim(pro_file,stringsAsFactors = FALSE)
    rnaseq <- read.delim(rna_file,stringsAsFactors = FALSE)
    
    #################### Annotate genes with chromosomes
    length(intersect(rownames(proteome), rownames(rnaseq)))
    # 11143 genes present in both proteome and rnaseq
    geneSyms <- union(rownames(proteome), rownames(rnaseq))
    # total 19487 unique genes in both set
    
    library(biomaRt)
    mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl",
                    host = "www.ensembl.org")
    gene_info <- getBM(attributes=c('hgnc_symbol', 'chromosome_name'), 
                       filters='hgnc_symbol', 
                       values=geneSyms, 
                       mart=mart)
    out_gene_file <- paste(out_dir,"/genes.tsv",sep = "")
    write.table(gene_info, out_gene_file, sep='\t', col.names=T, row.names=F)
    
    
    ## extract sex genes to partition the datasets
    gene_chromosome <- read.table(out_gene_file, header=T)
    sexgenes <- which(gene_chromosome$chromosome_name %in% c('X', 'Y'))
    sexgenes <- gene_chromosome$hgnc_symbol[sexgenes]
    
    
    cat("Impute data:", pro_file,"\n")
    
    proteome <- read.delim(pro_file,stringsAsFactors = FALSE)
    
    ##### partition proteome into sex genes and autosomal genes
    ## for sex genes, replace missing values as 0
    sexprot <- proteome[intersect(rownames(proteome), sexgenes), ]
    sexprot <- t(sexprot)
    sexprot[is.na(sexprot)] <- 0
    #write.table(sexprot, sprintf('outdir/for_2c/Testing_%d/sexprot.tsv', d), row.names=T, col.names=T, sep='\t')
    out_f1 <- paste(out_dir,"/sexprot.tsv",sep="")
    write.table(sexprot, out_f1, row.names=TRUE, col.names=TRUE, sep='\t')
    
    ## for autosomal genes, remove rows which > 50% of samples have missing value and impute missing value
    autopro <- proteome[setdiff(rownames(proteome), sexgenes), ]
    pmiss <- apply(autopro, 1, function(x) sum(is.na(x)))
    pmiss <- which(pmiss < 40)
    autopro <- autopro[pmiss, ]
    print(paste('dimension', dim(autopro)))
    
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    autopro <- t(autopro)
    autoimp <- missForest(autopro,parallelize = "variables")
    stopCluster(cl)
    #write.table(autoimp$ximp, sprintf('outdir/for_2c/Testing_%d/autoprot.tsv', d), row.names=T, col.names=T, sep='\t')
    out_f2 <- paste(out_dir,"/autoprot.tsv",sep = "")
    write.table(autoimp$ximp, out_f2, row.names=TRUE, col.names=TRUE, sep='\t')
    return(c(out_f1,out_f2,out_gene_file))
}

########## Obtain gene annotation with chromosome position (annotated in 1_preprocessing.R)
getSexGenes <- function(gene_file){
    #gene_chromosome <- read.table('kidney_genes.tsv', header=T)
    gene_chromosome <- read.table(gene_file, header=T)
    sexgenes <- which(gene_chromosome$chromosome_name %in% c('X', 'Y'))
    sexgenes <- gene_chromosome$hgnc_symbol[sexgenes]
    return(sexgenes)
}


########## Obtain clinical data
getClinical <- function(clifile){
    #clifile <- sprintf('../mislabeling/data/for_2c/Testing_%d/test_cli.tsv', d)
    #clinic  <- read.table(clifile, header=T, sep='\t', quote='')
    clinic  <- read.table(clifile, header=T, sep='\t', quote='')
    return(clinic)
}


########## Obtain RNAseq data
#### partition RNAseq into sex genes and autosomal genes
## remove rows with missing values in autosomal genes
## replace missing values as 0 in sex genes
getRNAseq <- function(rna_file,gene_file){
    rnaseq <- list()
    
    #rnafile <- sprintf('../mislabeling/data/for_2c/Testing_%d/test_rna.tsv', d)
    rnafile <- rna_file
    rnadata  <- read.table(rnafile, header=TRUE)
    colnames(rnadata) <- paste0('RNA_', 1:ncol(rnadata))
    sexgenes <- getSexGenes(gene_file)
    
    rna_atsm <- rnadata[setdiff(rownames(rnadata), sexgenes), ]
    rna_atsm <- na.exclude(rna_atsm)
    rna_atsm <- rna_atsm[rowSums(rna_atsm) != 0, ]
    rna_atsm <- t(rna_atsm)
    
    rna_sex  <- rnadata[intersect(rownames(rnadata), sexgenes), ]
    rna_sex[is.na(rna_sex)] <- 0
    rna_sex  <- rna_sex[rowSums(rna_sex) != 0, ]
    rna_sex  <- t(rna_sex)
    
    rnaseq$rna_atsm <- rna_atsm
    rnaseq$rna_sex  <- rna_sex
    
    return(rnaseq)
}


########## Obtain Imputed Proteomic data (imputed by 1_preprocessing.R)
getProteome <- function(s_file, ns_file){
    proteome <- list()
    
    #profile  <- sprintf('outdir/for_2c/Testing_%d/autoprot.tsv', d)
    profile <- ns_file
    pro_atsm <- read.table(profile, header=T)
    rownames(pro_atsm) <- paste0('PRO_', 1:nrow(pro_atsm))
    
    #profile  <- sprintf('outdir/for_2c/Testing_%d/sexprot.tsv', d)
    profile <- s_file
    pro_sex  <- read.table(profile, header=T)
    rownames(pro_sex) <- paste0('PRO_', 1:nrow(pro_sex))
    pro_sex  <- pro_sex[, colSums(pro_sex) != 0]
    
    proteome$pro_atsm <- as.matrix(pro_atsm)
    proteome$pro_sex  <- as.matrix(pro_sex)
    
    return(proteome)
}


########## Extract correlated genes (gene which correlation > 0.5)
extractCorGene <- function(rnamatrix, promatrix){
    cormatrix <- cor(as.matrix(rnamatrix), as.matrix(promatrix) )
    correlate <- diag(cormatrix)
    corgene <- names(which(correlate > 0.5))
    cat(sprintf('%d genes are highly correlated \n', length(corgene)))
    return(corgene)
}


########## Invert correlation matrix into probability matrix
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


########## Obtain the ranking order for each row (each RNA sample) of probability matrix
rankbyRow <- function(rankdist){
    malerank <- t(apply(rankdist, 1, function(x) rank(-x, ties.method='first')))
    colnames(malerank) <- colnames(rankdist)
    return(malerank)
}


########## Obtain the ranking order for each column (each PRO sample) of probability matrix
rankbyCol <- function(rankdist) {
    fmlerank <- apply(rankdist, 2, function(x) rank(-x, ties.method='first'))
    rownames(fmlerank) <- rownames(rankdist)
    return(fmlerank)
}


########## Perform stable matching algorithm on two matrices
#### malerank obtained by rankbyRow() function
#### fmlerank obtained by rankbyRow() function
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


########## Perform stable matching between RNA and PRO samples and return matching score
#getMatching <- function(probMatrix){
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


########## Obtain the weight of each training instances based on class distribution
getClassWeight <- function(labels){
    weight <- rep( 1, length(labels) )
    labelclass <- unique(labels)
    weight[labels == labelclass[1]] <- length(labels)/ 2 / sum(labels == labelclass[1])
    weight[labels == labelclass[2]] <- length(labels)/ 2 / sum(labels == labelclass[2])
    return(weight)
}


########## Training Regularized Logistic Regression model
trainGLM <- function(msiLabel, inputmtx, alpha){
    weight <- getClassWeight(msiLabel)
    if (sum(weight) != length(msiLabel)){
        cat('Error: sum(classweight) does not equal to length(msiLabel)! \n')
    }   # should be nrow(trainset)
    
    # perform cross validation of elasticnet to determine optimum lambda
    cv.glm <- cv.glmnet(as.matrix(inputmtx), msiLabel, family="binomial", weights=weight, alpha=alpha)
    (best_lambda <- cv.glm$lambda.1se)
    fit <- glmnet(as.matrix(inputmtx), msiLabel, family="binomial", weights=weight, alpha=alpha, lambda=best_lambda)
    
    return(fit)
}


########## Perform k-fold Cross Validation to predict clinical attributes
trainGLMcv <- function(msiLabel, inputmtx, alpha, k = 5){
    flds <- createFolds(msiLabel, k = k, list = TRUE, returnTrain = FALSE)
    predoutput <- data.frame(gender=msiLabel, gender_prob=0)
    
    for (f in 1:k) {
        testidx <- flds[[f]]
        numvar <- 0
        iter   <- 0
        while (numvar < 4 && iter < 50){
            fit1 <- trainGLM(msiLabel[-testidx], inputmtx[-testidx, ], 0.3)
            numvar1 <- sum(coef(fit1) > 0)
            if (numvar1 >= numvar) {
                fit <- fit1
                numvar <- numvar1
            }
            iter <- iter + 1
            cat('Fold', f, '- Iter', iter, '- Gender predictor model has', numvar, 'variable. \n')
        }
        predoutput$gender[testidx] <- predict(fit, as.matrix(inputmtx[testidx, ]), type='class')[, 1]
        predoutput$gender_prob[testidx] <- predict(fit, as.matrix(inputmtx[testidx, ]), type='response')[, 1]
    }
    return(predoutput)
}



########## Determine clinical swapping cases
determineClinicalSwapping <- function(traincli, cli_suspect){
    clinic_swap <- c()
    
    if (length(cli_suspect) <= 1){
        cat('No Clinical Swapping Cases Found! \n')
    } else  {
        cat('Clinical Label suspect swapping:', paste0(cli_suspect), '\n')
        subset <- traincli[cli_suspect, c('gender_prob', 'rgender_prob', 'pgender_prob', 'pred_gender')]
        
        clidist <- matrix(nrow=length(cli_suspect), ncol=length(cli_suspect))
        rownames(clidist) <- cli_suspect
        colnames(clidist) <- cli_suspect
        for (male in as.character(cli_suspect)) {
            for (fmle in as.character(cli_suspect)) {
                clidist[male, fmle] <- 1 - abs(subset[male, 'gender_prob'] - subset[fmle, 'pred_gender'])
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
    
    return(clinic_swap)
}


cmdargs <- commandArgs(TRUE)
pro_file <- cmdargs[1]
rna_file <- cmdargs[2]
sample_file <- cmdargs[3]
out_dir <- cmdargs[4]
#prefix <- cmdargs[5]

run_2c(pro_file, rna_file, sample_file, out_dir=out_dir)


