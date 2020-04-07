library(dplyr)

combine_methods=function(method1_folder, method2_folder, 
                         sample_annotation_file,
                         clinical_attributes = c("gender"), 
                         out_dir = "./", prefix = "test"){
  
  #sj_intermediate_file <- paste(method1_folder,"/intermediate.csv", sep = "")
  sj_intermediate_file <- paste(method1_folder,"/clinical_attributes_pred.tsv", sep = "")
  sj_correlation_file <- paste(method1_folder,"/correlation.csv", sep = "")
  modelA_result_file <- paste(method2_folder,"/test_ModelA_result.csv", sep = "")
  modelB_result_file <- paste(method2_folder,"/test_ModelB_result.csv", sep = "")
  
  cli_data_file <- sample_annotation_file
  cli_data <- read.delim(cli_data_file,stringsAsFactors = FALSE)
  #sjcli <- read.table(sprintf('../../../results/for_2b/SoonJye2/intermediate/2b_intermediate_%d.csv', d), sep=',', header=T)
  sjcli <- read.table(sj_intermediate_file, header=TRUE, sep = "\t", quote = FALSE)
  
  out_prefix <- prefix

  sentionpro <- read.table(modelA_result_file, sep=',', header=TRUE)
  sentionrna <- read.table(modelB_result_file, sep=',', header=TRUE)
  
  #corsample <- read.table(sprintf('../../../results/for_2b/SoonJye2/intermediate/2b_correlation_%d.csv',d,d), sep=',', header=T, row.names=1)
  corsample <- read.table(sj_correlation_file, sep=',', header=TRUE, row.names=1)
  corsample <- as.matrix(corsample)
  rankdist <- computeRankdist(corsample)
  
  malerank <- rankbyRow(rankdist)
  fmlerank <- rankbyCol(rankdist)
  
  matcher <- stableMarriage(malerank, fmlerank)
  #rnamatch <- as.numeric(sub('PRO_', '', matcher$malematch))
  #promatch <- as.numeric(sub('RNA_', '', matcher$fmlematch))
  rnamatch <- c()
  promatch <- c()
  for(i in 1:length(matcher$malematch)){
    rnamatch[i] <- which(rownames(corsample) == matcher$malematch[i])
  }
  
  for(i in 1:length(matcher$fmlematch)){
    promatch[i] <- which(rownames(corsample) == matcher$fmlematch[i])
  }
  
  total_samples <- nrow(sjcli)
  #nonmatch <- which(1:100 != rnamatch)
  nonmatch <- which(1:total_samples != rnamatch)
  length(nonmatch)
  
  
  ### get pairdist for shifting chain identification later
  #pairdist <- data.frame(rna=1:100, pro=rnamatch)
  pairdist <- data.frame(rna=1:total_samples, pro=rnamatch)
  pairdist$rnarank <- apply(pairdist, 1, function(x) malerank[x['rna'], x['pro']])
  pairdist$prorank <- apply(pairdist, 1, function(x) fmlerank[x['rna'], x['pro']])
  pairdist$distance <- pairdist$rnarank + pairdist$prorank
  pairdist$softmax <- apply(pairdist, 1, function(x) rankdist[x['rna'], x['pro']])
  
  
  ### Determine spurious match due to being left out (threshold is 5% of n)
  #if (sum(1:100 == rnamatch & pairdist$distance >= 4) > 0){
  if (sum(1:total_samples == rnamatch & pairdist$distance >= 4) > 0){
    #spumatch <- pairdist$rna[1:100 == rnamatch & pairdist$distance >= 4]
    spumatch <- pairdist$rna[1:total_samples == rnamatch & pairdist$distance >= 4]
    nonmatch <- c(nonmatch, spumatch)
    cat('Spurious match:', paste(spumatch), '\n')
    print(pairdist[spumatch,])
  }
  cat('\n Number of Non-matching =', length(nonmatch), '\n\n')
  
  
  #traincli <- sjcli[, c('sample', 'gender', 'msi')]
  traincli <- sjcli[, c('sample', clinical_attributes)]
  
  ######### Prediction probability #########
  cli_data_use <- cli_data %>% filter(sample %in% traincli$sample)
  cli_data_use <- cli_data_use[match(traincli$sample, cli_data_use$sample),]
  
  #traincli$gender_prob <- apply(traincli, 1, function(x) if (x['gender'] == 'Female') 0 else 1)
  # true probability
  # traincli$gender_prob <- apply(traincli, 1, function(x) if (x['gender'] == levels(traincli[,'gender'])[1]) 0 else 1)
  #traincli$msi_prob <- apply(traincli, 1, function(x) if (x['msi'] == 'MSI-High') 0 else 1)
  # true probability
  # traincli$msi_prob <- apply(traincli, 1, function(x) if (x['msi'] == levels(traincli[,'msi'])[1]) 0 else 1)
  
  # true probability
  cli_attr_prob_names_true <- paste(clinical_attributes,"_prob",sep = "")
  for(i in 1:length(clinical_attributes)){
    #cli_attr <- clinical_attributes[i]
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    traincli[,cli_attr_prob_names_true[i]] <- apply(traincli, 1, function(x) if (x[clinical_attributes[i]] == levels(traincli[,clinical_attributes[i]])[1]) 0 else 1)
  }
  
  
  cat('\n')
  ### Merge probability from three teams
  #traincli$rgender_prob <- (sjcli$GenderByRNA + lundcli$GenderByRNA + sentionrna$gender_prob) / 3
  #traincli$rmsi_prob <- (sjcli$MsiByRNA + lundcli$MsiByRNA + sentionrna$msi_prob) / 3
  #traincli$pgender_prob <- (sjcli$GenderByPRO + lundcli$GenderByPRO + sentionpro$gender_prob) / 3
  #traincli$pmsi_prob <- (sjcli$MsiByPRO + lundcli$MsiByPRO + sentionpro$msi_prob) / 3
  # RNA
  #traincli$rgender_prob <- (sjcli$GenderByRNA + sentionrna$gender_prob) / 2
  #traincli$rmsi_prob <- (sjcli$MsiByRNA + sentionrna$msi_prob) / 2
  # Protein
  #traincli$pgender_prob <- (sjcli$GenderByPRO + sentionpro$gender_prob) / 2
  #traincli$pmsi_prob <- (sjcli$MsiByPRO + sentionpro$msi_prob) / 2
  # combine RNA and Protein
  #traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob) / 2
  #traincli$pred_msi <- (traincli$rmsi_prob + traincli$pmsi_prob) / 2
  

  # prediction probability
  cli_attr_prob_names_pred_rna <- paste("r",clinical_attributes,"_prob",sep = "")
  cli_attr_prob_names_pred_pro <- paste("p",clinical_attributes,"_prob",sep = "")
  
  cli_attr_prob_names_pred_combine <- paste("pred_",clinical_attributes,sep = "")
  
  for(i in 1:length(clinical_attributes)){
    cat("clinical attributes: ",clinical_attributes[i],"\n")
    # RNA
    #traincli$rgender_prob <- (sjcli$GenderByRNA + sentionrna$gender_prob) / 2
    #traincli$rmsi_prob <- (sjcli$MsiByRNA + sentionrna$msi_prob) / 2
    traincli[,cli_attr_prob_names_pred_rna[i]] <- ( sjcli[, cli_attr_prob_names_pred_rna[i]] + sentionrna[, paste(clinical_attributes[i],"_prob",sep = "")] ) / 2.0
    
    # Protein
    #traincli$pgender_prob <- (sjcli$GenderByPRO + sentionpro$gender_prob) / 2
    #traincli$pmsi_prob <- (sjcli$MsiByPRO + sentionpro$msi_prob) / 2
    traincli[,cli_attr_prob_names_pred_pro[i]] <- ( sjcli[, cli_attr_prob_names_pred_pro[i]] + sentionpro[, paste(clinical_attributes[i],"_prob",sep = "")] ) / 2.0
    
    # combine RNA and Protein
    #traincli$pred_gender <- (traincli$rgender_prob + traincli$pgender_prob) / 2
    #traincli$pred_msi <- (traincli$rmsi_prob + traincli$pmsi_prob) / 2
    traincli[,cli_attr_prob_names_pred_combine[i]] <- ( traincli[,cli_attr_prob_names_pred_rna[i]] + traincli[,cli_attr_prob_names_pred_pro[i]]) / 2.0
  }
  
  
  
  ########### Label Correction
  #final_tab <- data.frame(sample=paste0('Testing_', 1:100), Clinical=1:100, RNAseq=1:100, Proteomics=1:100)
  final_tab <- data.frame(sample=cli_data_use$sample, Clinical=1:total_samples, RNAseq=1:total_samples, Proteomics=1:total_samples)
  
  
  cat('\n')
  ### Clinical Swapping
  # gender_prob is true probability and pred_gender is predicted probability
  
  ## only consider gender
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
    #cli_suspect <- which(abs(traincli$gender_prob - traincli$pred_gender) > 0.45)
    cli_suspect <- which(abs(traincli[,gender_prob] - traincli[,pred_gender]) > 0.45)
    cli_suspect <- setdiff(cli_suspect, nonmatch)
    traincli[cli_suspect, ]
  }else{
    cli_suspect <- c()
  }
  
  clinic_swap <- c()
  if (length(cli_suspect) <= 1){
    cat('No Clinical Swapping Cases Found! \n')
  } else  {
    cat('Suspected Clinical swap case:', paste(cli_suspect), '\n')
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
  
  clinic_swap <- union(clinic_swap, high_suspect)
  
  
  
  
  
  cat('\n')
  ### check swapping
  proswap <- c()
  rnaswap <- c()
  swapped <- c()
  for (r in nonmatch){
    if (r %in% swapped)    next
    p <- rnamatch[r]
    if (rnamatch[p] == r && rnamatch[r] != r) {
      if (malerank[r,p] + fmlerank[r,p] + malerank[p,r] + fmlerank[p,r] > 8){
        cat(sprintf('Spurious Pair: %d <--> %d, distance = %f \n', r, p, malerank[r,p] + fmlerank[r,p] + malerank[p,r] + fmlerank[p,r]))
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
  ### check duplication and shifting
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
    
    chains <- list()
    i <- 1
    for (start in lose_starts){
      chain <- c(start)
      cnext <- start
      while (!(cnext %in% lose_ends)) {
        cnext <- promatch[cnext]
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
      cat('Error! Not all suspected samples are found in chain. Circular shifting suspected! \n\n')
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
      } else{
        if (pro_shift < rna_shift) {          # supposedly distfront < distback
          cat('Proteome shift: ', chain[1], paste(chain[3:lenchain], collapse=' '), chain[length(chain)], sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: RNA_%d <-> PRO_%d = %.4f \t RNA_%d <-> PRO_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'Proteomics'] <- chain[3:length(chain)]
          proshift <- c(proshift, chain[2:(lenchain-1)])
          if (distfront > distback) {
            cat('Error: distance indicates RNA_', chain[1], ' duplication but prediction results indicate Proteome shifting', sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift), '\n')
            pdf(paste(out_dir,"/",out_prefix,"_pro_shifting.pdf",sep = ""),width = 4.5,height = 4.5)
            heatmap(rankdist[chain, chain], Rowv=NA, Colv=NA)
            dev.off()
          }
          
        } else if (rna_shift < pro_shift) {   # supposedly distfront > distback
          cat('RNA shift: ', chain[1], paste(chain[1:(length(chain)-2)], collapse = ' '), chain[length(chain)], sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift))
          cat(sprintf('Distance: RNA_%d <-> PRO_%d = %.4f \t RNA_%d <-> PRO_%d = %.4f \n', chain[2], chain[1], distfront, chain[length(chain)], chain[length(chain)-1], distback))
          final_tab[chain[-c(1,length(chain))], 'RNAseq'] <- chain[1:(length(chain)-2)]
          rnashift <- c(rnashift, chain[2:(lenchain-1)])
          if (distfront < distback){
            cat('Error: distance indicates PRO_', chain[lenchain], ' duplication but prediction results indicate RNAseq shifting', sprintf('(RNA: %.3f vs PRO: %.3f) \n', rna_shift, pro_shift), '\n')
            pdf(paste(out_dir,"/",out_prefix,"_rna_shifting.pdf",sep = ""),width = 4.5,height = 4.5)
            heatmap(rankdist[chain, chain], Rowv=NA, Colv=NA)
            dev.off()
          }
        }
      }
      cat('\n')
    }
  }
  
  cat('\n')
  out_file <- paste(out_dir,"/final_tab.tsv",sep = "")
  write.table(final_tab,out_file,col.names=TRUE, row.names=FALSE, sep=',')

}



