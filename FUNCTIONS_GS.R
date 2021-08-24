#A function that will filter a genotype matrix based on maf and missingness
#Calculates the proportion of missing data for every marker in a genotype matrix (mising data is NA)
calc_missrate <- function(gt_mat)
{
  col_func <- function(gt_col)
  {
    missrate <- sum(is.na(gt_col)) / length(gt_col)
    return(missrate)
  }

  missrate_vect <- apply(gt_mat, 2, col_func)

  return(missrate_vect)
}

# Calculates the minor allele frequency for every marker in a genotype matrix (coded as c(-1,0,1))
calc_maf_apply <- function(gt_mat, encoding = c(-1, 0, 1))
{
  col_func1 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == -1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
    allele2_ct <- (sum(gt_col == 1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)

    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }

  col_func2 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == 0, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
    allele2_ct <- (sum(gt_col == 2, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)

    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }

  if (all(encoding == c(-1, 0, 1)))
  {
    maf_vect <- apply(gt_mat, 2, col_func1)
  } else if (all(encoding == c(0, 1, 2)))
  {
    maf_vect <- apply(gt_mat, 2, col_func2)
  } else{
    print('Encoding not recognized, returning NULL')
    maf_vect <- NULL
  }

  return(maf_vect)
}

# This is a function that will split data into a list of k-folds
make_CV_sets <- function(list_length, k = 5)
{
  rand_values <- rnorm(list_length)
  k_quantiles <- quantile(rand_values, 0:k/k)
  k_assign <- cut(rand_values, k_quantiles, include.lowest = T, labels = F)

  cv_list <- list()
  for (i in 1:k)
  {
    fold_assignment <- k_assign != i
    cv_list[[i]] <- fold_assignment
  }
  return(cv_list)
}

test_all_models_BGLR_cv <- function(genotypes, phenotype, nIter = 5000, burnIn = 2000, folds = 5)
{
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = 5)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA

    # Calculate the GS model using BGLR
    ##Bayes A
    bayesA_ETA<-list(list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices])
    ##Bayes
    bayesB_ETA<-list(list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices])
    ##Bayes C
    bayesC_ETA<-list(list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices])
    ##Bayes RR
    BRR_ETA<-list(list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices])
    ##Bayes L
    BL_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices])
    ##Single Kernel RKHS
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices])
    ##Multi-Kernel RKHS
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices])
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }

  names(BGLR_acc_table) <- c("fold", "model", "r")

  result_plot <- ggplot(BGLR_acc_table, aes(x = model, y = r)) + geom_boxplot()

  results <- list(BGLR_acc_table, result_plot)
  names(results) <- c("BGLR_cv_acc_table", "BGLR_model_cv_boxplot")

  return(results)
}

test_all_models_BGLR_cv_mean <- function(genotypes, phenotype, nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = 5)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA

    # Calculate the GS model using BGLR
    ##Bayes A
    bayesA_ETA<-list(list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices])
    ##Bayes
    bayesB_ETA<-list(list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices])
    ##Bayes C
    bayesC_ETA<-list(list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices])
    ##Bayes RR
    BRR_ETA<-list(list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices])
    ##Bayes L
    BL_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices])
    ##Single Kernel RKHS
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices])
    ##Multi-Kernel RKHS
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices])
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:8])
  return(results)
}

rrblup_cv=function(x,m,p){
  library(rrBLUP)
  myCV_sets <- make_CV_sets(x,k = 5)
  pheno=data.frame(x)
  rrBLUP_results <- list()
  for (i in 1:length(myCV_sets)){
    train = which(myCV_sets[[i]])
    test = which(!myCV_sets[[i]])

    myGD_train <- as.matrix(m[train,])
    myGD_test <- as.matrix(m[test,])
    myY_train <- pheno[train,]
    myY_test <- pheno[test,]
    myPCA_train <- p[train,]
    myPCA_test <- p[test,]
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    results <- list(rrBLUP_model, predictions, acc)
    names(results) <- c("rrBLUP_model", "predictions", "acc")

    rrBLUP_results[[i]] <- results
  }
  rrBLUP_acc <- c(NULL)
  for (i in 1:length(myCV_sets)){
    rrBLUP_acc <- c(rrBLUP_acc, rrBLUP_results[[i]]$acc)
  }
  rrBLUP_acc}

rrblup_cv_mean=function(x,m,p){
  library(rrBLUP)
  myCV_sets <- make_CV_sets(x,k = 5)
  pheno=data.frame(x)
  rrBLUP_results <- list()
  for (i in 1:length(myCV_sets)){
    train = which(myCV_sets[[i]])
    test = which(!myCV_sets[[i]])

    myGD_train <- as.matrix(m[train,])
    myGD_test <- as.matrix(m[test,])
    myY_train <- pheno[train,]
    myY_test <- pheno[test,]
    myPCA_train <- p[train,]
    myPCA_test <- p[test,]
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    results <- list(rrBLUP_model, predictions, acc)
    names(results) <- c("rrBLUP_model", "predictions", "acc")

    rrBLUP_results[[i]] <- results
  }
  rrBLUP_acc <- c(NULL)
  for (i in 1:length(myCV_sets)){
    rrBLUP_acc <- c(rrBLUP_acc, rrBLUP_results[[i]]$acc)
  }
  mean(rrBLUP_acc)}

BGLR_cv_Ordinal <- function(genotypes, phenotype, nIter = 5000, burnIn = 2000, folds = 5)
{
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA

    # Calculate the GS model using BGLR
    ##Ordinal
    BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    BO_predictions <- predict(BO_results)
    BO_acc <- cor(as.numeric(phenotype[-fold_indices]), BO_predictions[-fold_indices])
    DF=BO_results$probs[-fold_indices,]
    probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
    probs=as.numeric(probs)
    BO_acc_cat=cor(as.numeric(phenotype[-fold_indices]), probs,use = "complete.obs")
    BGLR_acc_results[[i]] <- list(Acc=BO_acc,Acc_cat=BO_acc_cat)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("R_acc","Cat_acc")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  return(acc_fold)
}
PandG <- function(Pheno_Em_bl=NULL,myGD_EM=NULL,mytaxa_EM=NULL,myGM_EM=NULL){
  colnames(Pheno_Em_bl)[1]<-c("Genotype")
  rownames(Pheno_Em_bl) <- Pheno_Em_bl$Genotype

  Markers_Em_bl<-myGD_EM

  rownames(Markers_Em_bl) <- mytaxa_EM$V1
  colnames(Markers_Em_bl) <- myGM_EM$SNP


  Pheno_Em_bl <- Pheno_Em_bl[rownames(Pheno_Em_bl) %in% rownames(Markers_Em_bl),]
  Markers_Em_bl <- Markers_Em_bl[rownames(Markers_Em_bl) %in% rownames(Pheno_Em_bl),]

  Pheno_Em_bl <- Pheno_Em_bl[order(Pheno_Em_bl$Genotype),]
  Markers_Em_bl <- Markers_Em_bl[order(rownames(Markers_Em_bl)),]

  Pheno_Em_bl <- Pheno_Em_bl[order(Pheno_Em_bl$Genotype),]
  Markers_Em_bl <- Markers_Em_bl[order(rownames(Markers_Em_bl)),]

  myGM<-myGM_EM[myGM_EM$SNP %in% colnames(Markers_Em_bl),]
  myGD<-data.frame(rownames(Markers_Em_bl),Markers_Em_bl)
  colnames(myGD)[1]<-c("taxa")

  PCA_bl_em=prcomp(Markers_Em_bl)
  myPCA_bl_em=PCA_bl_em$x
  PGPCA=list(pheno=Pheno_Em_bl,geno=Markers_Em_bl,map=myGM,numeric=myGD,PC=myPCA_bl_em)
  return(PGPCA)
}

PGandCV <- function(Pheno_Em_bl=NULL,myGD_EM=NULL,mytaxa_EM=NULL,myGM_EM=NULL,CV=NULL){
  colnames(Pheno_Em_bl)[1]<-c("Genotype")
  rownames(Pheno_Em_bl) <- Pheno_Em_bl$Genotype

  Markers_Em_bl<-myGD_EM

  rownames(Markers_Em_bl) <- mytaxa_EM$V1
  colnames(Markers_Em_bl) <- myGM_EM$SNP

  colnames(CV)[1]<-c("Genotype")
  rownames(CV) <- CV$Genotype
  library(tidyr)
  library(Hmisc)
  for(i in 2:ncol(CV)){
    CV[,i]=impute(CV[,i])
    CV[,i]=as.numeric(CV[,i])
  }


  Pheno_Em_bl <- Pheno_Em_bl[rownames(Pheno_Em_bl) %in% rownames(Markers_Em_bl),]
  Markers_Em_bl <- Markers_Em_bl[rownames(Markers_Em_bl) %in% rownames(Pheno_Em_bl),]

  Pheno_Em_bl <- Pheno_Em_bl[order(Pheno_Em_bl$Genotype),]
  Markers_Em_bl <- Markers_Em_bl[order(rownames(Markers_Em_bl)),]

  CV <- CV[rownames(CV) %in% rownames(Pheno_Em_bl),]
  CV <- CV[order(rownames(CV)),]

  myGM<-myGM_EM[myGM_EM$SNP %in% colnames(Markers_Em_bl),]
  myGD<-data.frame(rownames(Markers_Em_bl),Markers_Em_bl)
  colnames(myGD)[1]<-c("taxa")

  PCA_bl_em=prcomp(Markers_Em_bl)
  myPCA_bl_em=PCA_bl_em$x
  PGPCA=list(pheno=Pheno_Em_bl,geno=Markers_Em_bl,map=myGM,numeric=myGD,PC=myPCA_bl_em,CV=CV)
  return(PGPCA)
}

#just need phenotypic vector and GWAS results
#Phenotypic vector
#Marker_Effects(Pheno=GBS_qam_adj19$pheno[,c(19)],GWAS=myGAPIT_BLINK_qam_adj_19)
Marker_Effects<-function(Pheno=NULL,GWAS=NULL,alpha=0.05,correction="Bonferonni",messages=TRUE){
  #Input your phenotypic vector
  myY_train <- Pheno
  if(messages==TRUE){print(paste0("Marker effects are being calculated and using the ",correction," correction."))}
  #Correction Method
  #Bonferonni Correction
  if(correction=="Bonferonni"){sig_markers <- which(GWAS$GWAS$P.value <= alpha/length(GWAS$GWAS$P.value))}
  #Regular P-value
  if(correction=="P-value"){sig_markers <- which(GWAS$GWAS$P.value <= alpha)}
  PC_mat <-GWAS$PCA[,-1]
  MAS_mat<-GWAS$GD[,-1]
  MAS_mat<-data.frame(MAS_mat[,sig_markers])
  if(length(sig_markers)==1){names(MAS_mat)<-GWAS$GWAS[sig_markers,]$SNP[[1]]}
  if(length(sig_markers)==0){print("Your GWAS is terrible and Halle thinks you suck.")}
  if(length(sig_markers)==0){stop("Your GWAS and correction method found zero significant markers.")}
  if(messages==TRUE){print(paste0("There are ",length(sig_markers)," significant marker(s)."))}
  GLM_data <- data.frame(myY_train,PC_mat,MAS_mat)
  names(GLM_data)[1] <- "Y"
  MAS_model <- lm(Y ~ ., data = GLM_data)
  Marker_Results=data.frame()
  cov_length=ncol(data.frame(myY_train))+ncol(PC_mat)
  for(i in cov_length+1:ncol(as.data.frame(MAS_mat)) ){
    #Null Model
    GLM_data_null <- data.frame(myY_train,PC_mat)
    names(GLM_data_null)[1] <- "Y"
    MAS_model_null <- lm(Y ~ ., data = GLM_data_null)
    R2null=summary(MAS_model_null)$r.squared

    #Model with 1 marker
    SNP_Name=names(GLM_data)[i]
    GLM_data_w1 <- data.frame(myY_train,PC_mat,GLM_data[,i])
    names(GLM_data_w1)[1] <- "Y"
    names(GLM_data_w1)[cov_length+1] <- SNP_Name
    MAS_model_w1 <- lm(Y ~ ., data = GLM_data_w1)
    R2w1=summary(MAS_model_w1)$r.squared
    R2w1m=R2w1-R2null #r2 for significant marker

    #Model without marker
    GLM_data_sm1 <- GLM_data[,-i]
    MAS_model_sm1 <- lm(Y ~ ., data = GLM_data_sm1)
    R2full=summary(MAS_model)$r.squared
    R2wo=summary(MAS_model_sm1)$r.squared
    R2wo1m=summary(MAS_model)$r.squared-summary(MAS_model_sm1)$r.squared #r2 for significant marker
    Full_Effects=MAS_model$coefficients[i]
    SNP_Effects=MAS_model_w1$coefficients[cov_length+1]

    results=data.frame(SNP_Name,R2null,R2w1,R2w1m,R2full,R2wo,R2wo1m,Full_Effects,SNP_Effects)
    names(results)=c("SNP_Name","R2.null.model","R2.add.single.marker","R2.of.single.marker","R2.full.model","R2.minus.single.marker","R2.diff","Effects.marker.full.model","Effects.single.marker")
    Marker_Results=rbind(Marker_Results,results)
  }
  if(messages==TRUE){print(paste0("Wow! Your signficant markers and covariates accounted for ",round(R2full*100,2),"% of the variation."))}
  return(Marker_Results)
}

test_all_models_mean <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices])
    ##Bayes
    gc()
    bayesB_ETA<-list(list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices])
    ##Bayes C
    gc()
    bayesC_ETA<-list(list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices])
    ##Bayes RR
    gc()
    BRR_ETA<-list(list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices])
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices])
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices])
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices])
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:10])
  return(results)
}

test_all_models_cv <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices])
    ##Bayes
    gc()
    bayesB_ETA<-list(list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices])
    ##Bayes C
    gc()
    bayesC_ETA<-list(list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices])
    ##Bayes RR
    gc()
    BRR_ETA<-list(list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices])
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices])
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices])
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices])
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  return(acc_fold)
}


test_all_models_cv_pc <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices])
    ##Bayes
    gc()
    bayesB_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices])
    ##Bayes C
    gc()
    bayesC_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices])
    ##Bayes RR
    gc()
    BRR_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices])
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices])
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(X=PCA,model="FIXED"),list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices])
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(X=PCA,model="FIXED"),
                 list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices])
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  return(acc_fold)
}

#Mean

test_all_models_mean <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices])
    ##Bayes
    gc()
    bayesB_ETA<-list(list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices])
    ##Bayes C
    gc()
    bayesC_ETA<-list(list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices])
    ##Bayes RR
    gc()
    BRR_ETA<-list(list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices])
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices])
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices])
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices])
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:10])
  return(results)
}

test_all_models_mean_recode <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)

    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)

    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]

    #PCA_bl_em=prcomp(myGD_train)
    #myPCA_train=PCA_bl_em$x

    #PCA_bl_em=prcomp(myGD_test)
    #myPCA_test=PCA_bl_em$x
    #myPCA_train=myPCA_train[,1:3]
    #myPCA_test=myPCA_test[,1:3]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes
    gc()
    bayesB_ETA<-list(list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes C
    gc()
    bayesC_ETA<-list(list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes RR
    gc()
    BRR_ETA<-list(list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices], use = "pairwise.complete")
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices], use = "pairwise.complete")
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(
                 list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices], use = "pairwise.complete")
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:10])
  return(results)
}

test_all_models_mean_pc <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices])
    ##Bayes
    gc()
    bayesB_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices])
    ##Bayes C
    gc()
    bayesC_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices])
    ##Bayes RR
    gc()
    BRR_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices])
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices])
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(X=PCA,model="FIXED"),list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices])
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(X=PCA,model="FIXED"),
                 list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices])
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:10])
  return(results)
}

test_all_models_mean_pc_recode <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)

    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)

    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]

    PCA_bl_train=prcomp(myGD_train)
    myPCA_train=PCA_bl_train$x

    PCA_bl_test=prcomp(myGD_test)
    myPCA_test=PCA_bl_test$x


    #PCA_bl_em=prcomp(myGD_train)
    #myPCA_train=PCA_bl_em$x

    #PCA_bl_em=prcomp(myGD_test)
    #myPCA_test=PCA_bl_em$x
    #myPCA_train=myPCA_train[,1:3]
    #myPCA_test=myPCA_test[,1:3]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesA"))
    bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn,saveAt='ba_')
    bayesA_predictions <- predict(bayesA_results)
    bayesA_acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes
    gc()
    bayesB_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesB"))
    bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn,saveAt='bb_')
    bayesB_predictions <- predict(bayesB_results)
    bayesB_acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes C
    gc()
    bayesC_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesC"))
    bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn,saveAt='bc_')
    bayesC_predictions <- predict(bayesC_results)
    bayesC_acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes RR
    gc()
    BRR_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BRR"))
    BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn,saveAt='brr_')
    BRR_predictions <- predict(BRR_results)
    BRR_acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices], use = "pairwise.complete")
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
    BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn,saveAt='bl_')
    BL_predictions <- predict(BL_results)
    BL_acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices], use = "pairwise.complete")
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(X=PCA,model="FIXED"),list(K=K,model='RKHS'))
    RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn,saveAt="RKHS_")
    RKHS_predictions <- predict(RKHS_results)
    RKHS_acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices], use = "pairwise.complete")
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(X=PCA,model="FIXED"),
                 list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn,saveAt="MK_")
    MK_predictions <- predict(MK_results)
    MK_acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices], use = "pairwise.complete")
    ##Ordinal
    #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    #BO_predictions <- predict(BO_results)
    #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

    BGLR_acc_results[[i]] <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[,2:10]
  return(results)
}



test_all_models_BLUP_pc_mean <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")

    BGLR_acc_results[[i]] <- list(acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}


test_all_models_BLUP_pc_mean_recode <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")

    BGLR_acc_results[[i]] <- list(acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}


test_all_models_BLUP_pc_mean_M <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    CV=impute(as.factor(CV))
    CV=as.numeric(CV)
    myCV_train <- CV[fold_indices]
    myCV_test <- CV[-fold_indices]
    fix_train <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test  <- as.matrix(cbind(myCV_test,myPCA_test))
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")

    BGLR_acc_results[[i]] <- list(acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}


test_all_models_BLUP_pc_mean_recode_M <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)

    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)

    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    CV=impute(as.factor(CV))
    CV=as.numeric(CV)
    myCV_train <- CV[fold_indices]
    myCV_test <- CV[-fold_indices]
    fix_train <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test  <- as.matrix(cbind(myCV_test,myPCA_test))
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")

    BGLR_acc_results[[i]] <- list(acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}


test_all_models_BLUP_cv <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- myPCA_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    ##GBLUP
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    gBLUP_model <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- gBLUP_model$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    g_acc <- cor(predictions, myY_test, use = "pairwise.complete")

    BGLR_acc_results[[i]] <- list(acc,g_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  return(acc_fold)
}

test_all_models_vs_pc_mean <- function(train_genotypes, train_phenotype,train_PCA=NULL,test_genotypes, test_phenotype,test_PCA=NULL,nIter = 5000, burnIn = 2000)
{
  library(BGLR)
  library(tidyr)
  phenotype <- c(train_phenotype,test_phenotype)
  pheno_train <- c(train=train_phenotype,test=test_phenotype)
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)

  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_test <- as.matrix(test_genotypes)
  myY_train <- train_phenotype
  myY_test <- test_phenotype

  myPCA_train=train_PCA

  myPCA_test=test_PCA


  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = myPCA_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- myPCA_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  ##GBLUP
  gc()
  K <- tcrossprod(as.matrix(genotypes))
  #GBlup with fixed effects
  gBLUP_model <- mixed.solve(y = pheno_train,
                             K = K,
                             X=PCA)
  pred_effects <- as.matrix(gBLUP_model$u)[-c(1:length(train_phenotype)),]
  fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
  # Calculate the GS model using BGLR
  ##Bayes A
  gc()
  bayesA_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesA"))
  bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
  bayesA_predictions <- predict(bayesA_results)
  bayesA_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesA_predictions[-c(1:length(train_phenotype))])
  ##Bayes
  gc()
  bayesB_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesB"))
  bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
  bayesB_predictions <- predict(bayesB_results)
  bayesB_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesB_predictions[-c(1:length(train_phenotype))])
  ##Bayes C
  gc()
  bayesC_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesC"))
  bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
  bayesC_predictions <- predict(bayesC_results)
  bayesC_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesC_predictions[-c(1:length(train_phenotype))])
  ##Bayes RR
  gc()
  BRR_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BRR"))
  BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
  BRR_predictions <- predict(BRR_results)
  BRR_acc <- cor(phenotype[-c(1:length(train_phenotype))], BRR_predictions[-c(1:length(train_phenotype))])
  ##Bayes L
  gc()
  BL_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
  BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
  BL_predictions <- predict(BL_results)
  BL_acc <- cor(phenotype[-c(1:length(train_phenotype))], BL_predictions[-c(1:length(train_phenotype))])
  ##Single Kernel RKHS
  gc()
  X=as.matrix(genotypes)
  p<-ncol(X)
  X<-scale(X,center=TRUE,scale=TRUE)
  D<-(as.matrix(dist(X,method='euclidean'))^2)/p
  h<-0.5
  K<-exp(-h*D)
  RKHS_ETA<-list(list(X=PCA,model="FIXED"),list(K=K,model='RKHS'))
  RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
  RKHS_predictions <- predict(RKHS_results)
  RKHS_acc <- cor(phenotype[-c(1:length(train_phenotype))], RKHS_predictions[-c(1:length(train_phenotype))])
  ##Multi-Kernel RKHS
  gc()
  X=as.matrix(genotypes)
  X<-scale(X,center=TRUE,scale=TRUE)
  D<-(as.matrix(dist(X,method='euclidean'))^2)/p
  h<-0.5*c(1/5,1,5)
  MK_ETA<-list(list(X=PCA,model="FIXED"),
               list(K=exp(-h[1]*D),model='RKHS'),
               list(K=exp(-h[2]*D),model='RKHS'),
               list(K=exp(-h[3]*D),model='RKHS'))
  MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
  MK_predictions <- predict(MK_results)
  MK_acc <- cor(phenotype[-c(1:length(train_phenotype))], MK_predictions[-c(1:length(train_phenotype))])
  ##Ordinal
  #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
  #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
  #BO_predictions <- predict(BO_results)
  #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

  BGLR_acc_results <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)

  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(BGLR_acc_results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:10])
  return(results)
}


test_all_models_vs_pc_mean_recode <- function(train_genotypes, train_phenotype,train_PCA=NULL,test_genotypes, test_phenotype,test_PCA=NULL,nIter = 5000, burnIn = 2000)
{
  library(BGLR)
  library(tidyr)
  phenotype <- c(train_phenotype,test_phenotype)
  pheno_train <- c(train=train_phenotype,test=test_phenotype)
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)
  # Calculate the GS model using rrBLUP

  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  myY_train <- train_phenotype
  myY_test <- test_phenotype

  myPCA_train=train_PCA
  myPCA_test=test_PCA


  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = myPCA_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- myPCA_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  ##GBLUP
  gc()
  K <- tcrossprod(as.matrix(genotypes))
  #GBlup with fixed effects
  gBLUP_model <- mixed.solve(y = pheno_train,
                             K = K,
                             X=PCA)
  pred_effects <- as.matrix(gBLUP_model$u)[-c(1:length(train_phenotype)),]
  fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
  # Calculate the GS model using BGLR
  ##Bayes A
  gc()
  bayesA_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesA"))
  bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
  bayesA_predictions <- predict(bayesA_results)
  bayesA_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesA_predictions[-c(1:length(train_phenotype))])
  ##Bayes
  gc()
  bayesB_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesB"))
  bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
  bayesB_predictions <- predict(bayesB_results)
  bayesB_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesB_predictions[-c(1:length(train_phenotype))])
  ##Bayes C
  gc()
  bayesC_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesC"))
  bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
  bayesC_predictions <- predict(bayesC_results)
  bayesC_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesC_predictions[-c(1:length(train_phenotype))])
  ##Bayes RR
  gc()
  BRR_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BRR"))
  BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
  BRR_predictions <- predict(BRR_results)
  BRR_acc <- cor(phenotype[-c(1:length(train_phenotype))], BRR_predictions[-c(1:length(train_phenotype))])
  ##Bayes L
  gc()
  BL_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
  BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
  BL_predictions <- predict(BL_results)
  BL_acc <- cor(phenotype[-c(1:length(train_phenotype))], BL_predictions[-c(1:length(train_phenotype))])
  ##Single Kernel RKHS
  gc()
  X=as.matrix(genotypes)
  p<-ncol(X)
  X<-scale(X,center=TRUE,scale=TRUE)
  D<-(as.matrix(dist(X,method='euclidean'))^2)/p
  h<-0.5
  K<-exp(-h*D)
  RKHS_ETA<-list(list(X=PCA,model="FIXED"),list(K=K,model='RKHS'))
  RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
  RKHS_predictions <- predict(RKHS_results)
  RKHS_acc <- cor(phenotype[-c(1:length(train_phenotype))], RKHS_predictions[-c(1:length(train_phenotype))])
  ##Multi-Kernel RKHS
  gc()
  X=as.matrix(genotypes)
  X<-scale(X,center=TRUE,scale=TRUE)
  D<-(as.matrix(dist(X,method='euclidean'))^2)/p
  h<-0.5*c(1/5,1,5)
  MK_ETA<-list(list(X=PCA,model="FIXED"),
               list(K=exp(-h[1]*D),model='RKHS'),
               list(K=exp(-h[2]*D),model='RKHS'),
               list(K=exp(-h[3]*D),model='RKHS'))
  MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
  MK_predictions <- predict(MK_results)
  MK_acc <- cor(phenotype[-c(1:length(train_phenotype))], MK_predictions[-c(1:length(train_phenotype))])
  ##Ordinal
  #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
  #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
  #BO_predictions <- predict(BO_results)
  #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

  BGLR_acc_results <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)

  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(BGLR_acc_results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:10])
  return(results)
}


test_all_models_vs_mean_recode <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,nIter = 5000, burnIn = 2000)
{
  library(BGLR)
  library(tidyr)
  phenotype <- c(train_phenotype,test_phenotype)
  pheno_train <- c(train=train_phenotype,test=test_phenotype)
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  # Calculate the GS model using rrBLUP

  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  myY_train <- train_phenotype
  myY_test <- test_phenotype



  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  ##GBLUP
  gc()
  K <- tcrossprod(as.matrix(genotypes))
  #GBlup with fixed effects
  gBLUP_model <- mixed.solve(y = pheno_train,
                             K = K)
  pred_effects <- as.matrix(gBLUP_model$u)[-c(1:length(train_phenotype)),]
  fix_effects <-  gBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  g_acc <- cor(predictions, myY_test, use = "pairwise.complete")
  # Calculate the GS model using BGLR
  ##Bayes A
  gc()
  bayesA_ETA<-list(list(X=as.matrix(genotypes),model="BayesA"))
  bayesA_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn)
  bayesA_predictions <- predict(bayesA_results)
  bayesA_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesA_predictions[-c(1:length(train_phenotype))])
  ##Bayes
  gc()
  bayesB_ETA<-list(list(X=as.matrix(genotypes),model="BayesB"))
  bayesB_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn)
  bayesB_predictions <- predict(bayesB_results)
  bayesB_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesB_predictions[-c(1:length(train_phenotype))])
  ##Bayes C
  gc()
  bayesC_ETA<-list(list(X=as.matrix(genotypes),model="BayesC"))
  bayesC_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn)
  bayesC_predictions <- predict(bayesC_results)
  bayesC_acc <- cor(phenotype[-c(1:length(train_phenotype))], bayesC_predictions[-c(1:length(train_phenotype))])
  ##Bayes RR
  gc()
  BRR_ETA<-list(list(X=as.matrix(genotypes),model="BRR"))
  BRR_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn)
  BRR_predictions <- predict(BRR_results)
  BRR_acc <- cor(phenotype[-c(1:length(train_phenotype))], BRR_predictions[-c(1:length(train_phenotype))])
  ##Bayes L
  gc()
  BL_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
  BL_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn)
  BL_predictions <- predict(BL_results)
  BL_acc <- cor(phenotype[-c(1:length(train_phenotype))], BL_predictions[-c(1:length(train_phenotype))])
  ##Single Kernel RKHS
  gc()
  X=as.matrix(genotypes)
  p<-ncol(X)
  X<-scale(X,center=TRUE,scale=TRUE)
  D<-(as.matrix(dist(X,method='euclidean'))^2)/p
  h<-0.5
  K<-exp(-h*D)
  RKHS_ETA<-list(list(K=K,model='RKHS'))
  RKHS_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn)
  RKHS_predictions <- predict(RKHS_results)
  RKHS_acc <- cor(phenotype[-c(1:length(train_phenotype))], RKHS_predictions[-c(1:length(train_phenotype))])
  ##Multi-Kernel RKHS
  gc()
  X=as.matrix(genotypes)
  X<-scale(X,center=TRUE,scale=TRUE)
  D<-(as.matrix(dist(X,method='euclidean'))^2)/p
  h<-0.5*c(1/5,1,5)
  MK_ETA<-list(
               list(K=exp(-h[1]*D),model='RKHS'),
               list(K=exp(-h[2]*D),model='RKHS'),
               list(K=exp(-h[3]*D),model='RKHS'))
  MK_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn)
  MK_predictions <- predict(MK_results)
  MK_acc <- cor(phenotype[-c(1:length(train_phenotype))], MK_predictions[-c(1:length(train_phenotype))])
  ##Ordinal
  #BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
  #BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
  #BO_predictions <- predict(BO_results)
  #BO_acc <- cor(phenotype[-fold_indices], BO_predictions[-fold_indices])

  BGLR_acc_results <- list(acc,g_acc,bayesA_acc, bayesB_acc, bayesC_acc, BRR_acc, BL_acc, RKHS_acc, MK_acc)

  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP","BayesA", "BayesB", "BayesC", "BRR", "BL", "RKHS", "MK")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(BGLR_acc_results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:10])
  return(results)
}

test_BLUPs_vs_pc_mean <- function(train_genotypes, train_phenotype,train_PCA=NULL,test_genotypes, test_phenotype,test_PCA=NULL,nIter = 5000, burnIn = 2000)
{
  library(BGLR)
  library(tidyr)
  pheno_train <- c(train=train_phenotype,test=test_phenotype)
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)
  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_test <- as.matrix(test_genotypes)
  myY_train <- train_phenotype
  myY_test <- test_phenotype
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA

  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = myPCA_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- myPCA_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  ##GBLUP
  gc()
  K <- tcrossprod(as.matrix(genotypes))
  #GBlup with fixed effects
  gBLUP_model <- mixed.solve(y = pheno_train,
                             K = K,
                             X=PCA)
  pred_effects <- as.matrix(gBLUP_model$u)[-c(1:length(train_phenotype)),]
  fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  g_acc <- cor(predictions, myY_test, use = "pairwise.complete")

  BGLR_acc_results <- list(acc,g_acc)

  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(BGLR_acc_results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("VS", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:3])
  return(results)
}

test_BLUPs_vs_pc_mean_recode <- function(train_genotypes, train_phenotype,train_PCA=NULL,test_genotypes, test_phenotype,test_PCA=NULL,nIter = 5000, burnIn = 2000)
{
  library(BGLR)
  library(tidyr)
  train_genotypes = GBS_qam_adj18$geno
  train_phenotype = GBS_qam_adj18$pheno[,18]
  train_PCA=GBS_qam_adj18$PC[,1:3]
  test_genotypes=GBS_qam_adj19$geno
  test_phenotype = GBS_qam_adj19$pheno[,19]
  test_PCA=GBS_qam_adj19$PC[,1:3]

  pheno_train <- c(train=train_phenotype,test=test_phenotype)
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)
  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  myY_train <- train_phenotype
  myY_test <- test_phenotype
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA

  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = myPCA_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- myPCA_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  ##GBLUP
  gc()
  K <- tcrossprod(as.matrix(genotypes))
  #GBlup with fixed effects
  gBLUP_model <- mixed.solve(y = pheno_train,
                             K = K,
                             X=PCA)
  pred_effects <- as.matrix(gBLUP_model$u)[-c(1:length(train_phenotype)),]
  fix_effects <- as.matrix(myPCA_test)  %*% gBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  g_acc <- cor(predictions, myY_test, use = "pairwise.complete")

  BGLR_acc_results <- list(acc,g_acc)

  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP","gBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(BGLR_acc_results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("VS", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:3])
  return(results)
}

BLUP_A_pc_mean <- function(genotypes, phenotype,GM, folds = 5)
{
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])



    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices,2] <- NA
    pheno_test=phenotype[-fold_indices,]
    # Calculate the GS model using rrBLUP
    names(pheno_train)=c("Taxa", "Phenotype")
    names(pheno_test)=c("Taxa", "Phenotype")
    names(genotypes)[1]=c("Taxa")

    gc()

    myGAPIT6=GAPIT(
      Y=pheno_train,
      GD=genotypes,
      GM=GM,
      PCA.total=3,
      model="cBLUP",
      group.from=100,
      group.to=1000,
      group.by=20,
      file.output=F)
    gapit6=merge(pheno_test,myGAPIT6$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
    c_acc=cor(gapit6[,2],gapit6[,5])
    gc()
    myGAPIT7=GAPIT(
      Y=pheno_train,
      GD=genotypes,
      GM=GM,
      PCA.total=3,
      model="sBLUP",
      file.output=F)
    gapit7=merge(pheno_test,myGAPIT7$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
    s_acc=cor(gapit7[,2],gapit7[,5])

    BGLR_acc_results[[i]] <- list(c_acc,s_acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("cBLUP","sBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2:3]
  return(results)
}

BLUP_A_vs_pc_mean <- function(train_genotypes, train_phenotype,train_GM=NULL,test_genotypes, test_phenotype,test_GM=NULL, burnIn = 2000)
{
  library(tidyr)

  names(train_phenotype)=c("Taxa", "Phenotype")
  names(test_phenotype)=c("Taxa", "Phenotype")

  pheno_test=test_phenotype
  pheno_test[,2]<-NA
  pheno_train <- rbind(train_phenotype,pheno_test)

  genotypes<-rbind(train_genotypes,test_genotypes)

  names(pheno_train)=c("Taxa", "Phenotype")
  names(test_phenotype)=c("Taxa", "Phenotype")
  names(genotypes)[1]=c("Taxa")
  gc()

  myGAPIT6=GAPIT(
    Y=pheno_train,
    GD=genotypes,
    GM=train_GM,
    PCA.total=3,
    model="cBLUP",
    group.from=100,
    group.to=1000,
    group.by=20,
    file.output=F)
  gapit6=merge(test_phenotype,myGAPIT6$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
  c_acc=cor(gapit6[,2],gapit6[,5])
  gc()
  myGAPIT7=GAPIT(
    Y=pheno_train,
    GD=genotypes,
    GM=train_GM,
    PCA.total=3,
    model="sBLUP",
    file.output=F)
  gapit7=merge(test_phenotype,myGAPIT7$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
  s_acc=cor(gapit7[,2],gapit7[,5])


  BGLR_acc_results <- list(c_acc,s_acc)

  #, GBLUP_acc, "GBLUP"
  model_vect <- c("cBLUP","sBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(BGLR_acc_results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2:3]
  return(results)
}


Caret_Models_Mean_M <- function(genotypes, phenotype,CV=NULL,model="rf", folds = 5,markers=3000){
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myY_train=droplevels(myY_train)
    samp=sample(1:length(genotypes), markers)
    m_samp=genotypes[,samp]
    myGD_train <- m_samp[fold_indices,]
    myGD_test <- m_samp[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    if(length(mono_indices)==0){
      myGD_train = myGD_train
      myGD_test = myGD_test

    }else{

      myGD_train = myGD_train[,-mono_indices]
      myGD_test = myGD_test[,-mono_indices]
    }
    mydata <- data.frame(myY_train, myGD_train)
    mydata$myY_train=as.factor(mydata$myY_train)
    myY_test=as.factor(myY_test)
    colnames(mydata)[1]<-c("Y")
    gc()
    svmFit1 <- train(Y ~ ., data = mydata,
                     method = model,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 10),
                     tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myGD_test)
    acc=postResample(pred = svm.linear_pred1, obs = myY_test)
    #cm=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
    #acc=acc[1]
    svm_results[[i]] <- list(acc[1],acc[2])
  }
  model_vect <- c('r2','kappa')
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2:3]
  return(results)
}

Caret_Models_Mean_M_CV <- function(genotypes, phenotype,CV=NULL,model="rf", folds = 5,markers=3000)
{
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    myY_train=droplevels(myY_train)
    samp=sample(1:length(genotypes), markers)
    m_samp=genotypes[,samp]
    myGD_train <- m_samp[fold_indices,]
    myGD_test <- m_samp[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    if(length(mono_indices)==0){
      myGD_train = myGD_train
      myGD_test = myGD_test

    }else{

      myGD_train = myGD_train[,-mono_indices]
      myGD_test = myGD_test[,-mono_indices]
    }

    CV=impute(as.factor(CV))
    CV=as.numeric(CV)
    myCV_train <- CV[fold_indices]
    myCV_test <- CV[-fold_indices]
    fix_train <- cbind(myCV_train,myPCA_train)
    fix_test  <- cbind(myCV_test,myPCA_test)
    mydata <- data.frame(myY_train,myFix_train)
    mydata$myY_train=as.factor(mydata$myY_train)
    myY_test=as.factor(myY_test)
    colnames(mydata)[1]<-c("Y")
    gc()
    svmFit1 <- train(Y ~ ., data = mydata,
                     method = model,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 10),
                     tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myFix_test)
    acc=postResample(pred = svm.linear_pred1, obs = myY_test)
    #cm=confusionMatrix(data = svm.linear_pred1, reference = myY_test)
    acc=acc[1]
    results <- list(acc)
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

Caret_Models_Regression <- function(genotypes, phenotype,model="rf", folds = 5,markers=3000)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    samp=sample(1:length(genotypes), markers)
    m_samp=genotypes[,samp]
    myGD_train <- m_samp[fold_indices,]
    myGD_test <- m_samp[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    myGD_train = myGD_train[,-mono_indices]
    myGD_test = myGD_test[,-mono_indices]
    mydata <- data.frame(myY_train, myGD_train)
    colnames(mydata)[1]<-c("Y")
    predVars <- names(mydata)[-1]
    gc()
    svmFit1 <- train(mydata[,"Y"],
                     mydata[,predVars],
                     method = model)
    svm.linear_pred1 <- predict(svmFit1, myGD_test)
    acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

    results <- list(acc)
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  return(acc_fold)
}

Caret_Models_Regression_Mean_PC <- function(PCA, phenotype,model="svmRadial", folds = 5)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)# for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]


    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]


    mydata <- data.frame(myY_train, myPCA_train)
    colnames(mydata)[1]<-c("Y")

    gc()
    svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                     data = mydata,
                     method = model,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 10),
                     tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myPCA_test)
    acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

    results <- list(acc)
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

Caret_Models_Regression_Mean_M <- function(genotypes, phenotype,model="svmRadial", folds = 5,markers=3000)
{
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    samp=sample(1:length(genotypes), markers)
    m_samp=genotypes[,samp]
    myGD_train <- m_samp[fold_indices,]
    myGD_test <- m_samp[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    if(length(mono_indices)==0){
      myGD_train = myGD_train
      myGD_test = myGD_test
    }else{
    myGD_train = myGD_train[,-mono_indices]
    myGD_test = myGD_test[,-mono_indices]}
    mydata <- data.frame(myY_train, myGD_train)
    colnames(mydata)[1]<-c("Y")
    svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                     data = mydata,
                     method = model,
                     trControl = trainControl(method = "cv", number = 10),
                     tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myGD_test)
    acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

    results <- acc
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

RF_Mean_PC <- function(PCA, phenotype,model="rf", folds = 5)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)
  library(randomForest)# for classification and regression training
  library(kernlab)
  library(ranger)
  library(h2o)
  library(e1071)# for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]


    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]

    mydata <- data.frame(myY_train, myPCA_train)
    colnames(mydata)[1]<-c("Y")
    gc()
    svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                     data = mydata,
                     method = model,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 10),
                     tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myPCA_test)
    acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

    results <- list(acc)
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

RF_Mean_M <- function(genotypes, phenotype,model="rf", folds = 5,markers=3000){
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)
  library(randomForest)
  library(ranger)
  library(h2o)
  library(e1071)  # for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    samp=sample(1:length(genotypes), markers)
    m_samp=genotypes[,samp]
    myGD_train <- m_samp[fold_indices,]
    myGD_test <- m_samp[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    if(length(mono_indices)==0){
      myGD_train = myGD_train
      myGD_test = myGD_test
    }else{
      myGD_train = myGD_train[,-mono_indices]
      myGD_test = myGD_test[,-mono_indices]}
    mydata <- data.frame(myY_train, myGD_train)
    colnames(mydata)[1]<-c("Y")
    gc()
    svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                     data = mydata,
                     method = model,
                     preProcess = c("center", "scale"),
                     trControl = trainControl(method = "cv", number = 10),
                     tuneLength = 10)
    svm.linear_pred1 <- predict(svmFit1, myGD_test)
    acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

    results <- list(acc)
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}


NN_Mean_M <- function(genotypes, phenotype,model="neuralnet", folds = 5,markers=3000){
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)
  library(randomForest)
  library(ranger)
  library(h2o)
  library(e1071)
  library(neuralnet)
  # for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]
    samp=sample(1:length(genotypes), markers)
    m_samp=genotypes[,samp]
    myGD_train <- m_samp[fold_indices,]
    myGD_test <- m_samp[-fold_indices,]

    maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
    mono_indices <- which(maf == 0)
    #taxa <- row.names(myGD)
    if(length(mono_indices)==0){
      myGD_train = myGD_train
      myGD_test = myGD_test
    }else{
      myGD_train = myGD_train[,-mono_indices]
      myGD_test = myGD_test[,-mono_indices]}

    gc()
    tooHigh <- findCorrelation(cor(myGD_train), cutoff = .90)
    myGD_train = myGD_train[,-tooHigh]
    myGD_test = myGD_test[,-tooHigh]
    mydata <- data.frame(myY_train, myGD_train)
    colnames(mydata)[1]<-c("Y")
    ## Create a specific candidate set of models to evaluate:
    svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                     data = mydata,
                     method = "neuralnet",
                     trControl = trainControl(method = "cv", number = 10),
                     ## Automatically standardize data prior to modeling
                     ## and prediction
                     preProc = c("center", "scale"),
                     linout = TRUE,
                     trace = FALSE,
                     MaxNWts = 10 * (ncol(myGD_train) + 1) + 10 + 1,
                     maxit = 500)
    svm.linear_pred1 <- predict(svmFit1, myGD_test)
    acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

    results <- list(acc)
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

NN_Mean_PC <- function(PCA, phenotype,model="avNNet", folds = 5)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)
  library(randomForest)# for classification and regression training
  library(kernlab)
  library(ranger)
  library(h2o)
  library(e1071)# for fitting SVMs
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)
  svm_results <- list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]

    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])

    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]

    gc()
    tooHigh <- findCorrelation(cor(myPCA_train), cutoff = .90)
    myPCA_train <- myPCA_train[-tooHigh,]
    myPCA_test <- myPCA_test[-tooHigh,]
    mydata <- data.frame(myY_train, myPCA_train)
    colnames(mydata)[1]<-c("Y")
    ## Create a specific candidate set of models to evaluate:
    nnetGrid <- expand.grid(decay = c(0, 0.01, .1),
                            size = c(1:10),
                            ## The next option is to use bagging (see the
                            ## next chapter) instead of different random
                            ## seeds.
                            bag = FALSE)
    svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                     data = mydata,
                     method = model,
                     tuneGrid = nnetGrid,
                     trControl = trainControl(method = "cv", number = 10),
                     ## Automatically standardize data prior to modeling
                     ## and prediction
                     preProc = c("center", "scale"),
                     linout = TRUE,
                     trace = FALSE,
                     MaxNWts = 10 * (ncol(myPCA_train) + 1) + 10 + 1,
                     maxit = 500)
    svm.linear_pred1 <- predict(svmFit1, myPCA_test)
    acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

    results <- list(acc)
    names(results) <- c("acc")
    svm_results[[i]] <- results
  }
  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

SVM_Mean_VS_PC <- function(train_PCA, train_phenotype,test_PCA,test_phenotype,model="svmRadial", folds = 5)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  myY_train <- train_phenotype
  myY_test <- test_phenotype


  if(ncol(test_PCA)<ncol(train_PCA)){
    myPCA_train <- train_PCA
    myPCA_test <- test_PCA

  }else{
    myPCA_train <- train_PCA[,c(1:10)]
    myPCA_test <- test_PCA[,c(1:10)]
  }


  mydata <- data.frame(myY_train, myPCA_train)
  colnames(mydata)[1]<-c("Y")
  gc()
  svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                   data = mydata,
                   method = model,
                   preProcess = c("center", "scale"),
                   trControl = trainControl(method = "cv", number = 10),
                   tuneLength =10)

  svm.linear_pred1 <- predict(svmFit1, newdata =myPCA_test)
  acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")


  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  SVM_acc_table <- rbind(SVM_acc_table, results_long)

  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

SVM_Mean_VS_M <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,model="svmRadial", folds = 5,markers=3000)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)  # for fitting SVMs
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  myY_train <- train_phenotype
  myY_test <- test_phenotype
  samp=sample(1:length(train_genotypes), markers)

  myGD_train <- train_genotypes[,samp]
  myGD_test <- test_genotypes[,samp]

  maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
  mono_indices <- which(maf == 0)

  if(length(mono_indices)==0){
    myGD_train = myGD_train
    myGD_test = myGD_test
  }else{
    myGD_train = myGD_train[,-mono_indices]
    myGD_test = myGD_test[,-mono_indices]}

  mydata <- data.frame(myY_train, myGD_train)
  colnames(mydata)[1]<-c("Y")
  gc()
  svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                   data = mydata,
                   method = model,
                   preProcess = c("center", "scale"),
                   trControl = trainControl(method = "cv", number = 10),
                   tuneLength = 10)
  svm.linear_pred1 <- predict(svmFit1, myGD_test)
  acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")

  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  SVM_acc_table <- rbind(SVM_acc_table, results_long)

  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

RF_Mean_VS_PC <- function(train_PCA, train_phenotype,test_PCA,test_phenotype,model="nnet", folds = 5)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)
  library(randomForest)# for classification and regression training
  library(kernlab)
  library(ranger)
  library(h2o)
  library(e1071)# for fitting SVMs
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  myY_train <- train_phenotype
  myY_test <- test_phenotype

  if(ncol(test_PCA)<ncol(train_PCA)){
    myPCA_train <- train_PCA[,c(1:ncol(test_PCA))]
    myPCA_test <- test_PCA

  }else{
    myPCA_train <- train_PCA
    myPCA_test <- test_PCA[,c(1:ncol(train_PCA))]
  }

  mydata <- data.frame(myY_train, myPCA_train)
  colnames(mydata)[1]<-c("Y")
  gc()
  svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                   data = mydata,
                   method = model,
                   preProcess = c("center", "scale"),
                   trControl = trainControl(method = "cv", number = 10),
                   tuneLength = 10)
  svm.linear_pred1 <- predict(svmFit1, myPCA_test)
  acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")

  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  SVM_acc_table <- rbind(SVM_acc_table, results_long)
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

RF_Mean_VS_M <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,model="nnet", folds = 5,markers=3000){
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)
  library(randomForest)
  library(ranger)
  library(h2o)
  library(e1071)  # for fitting SVMs
  # Make the CV list
  myY_train <- train_phenotype
  myY_test <- test_phenotype
  samp=sample(1:length(train_genotypes), markers)

  myGD_train <- train_genotypes[,samp]
  myGD_test <- test_genotypes[,samp]

  maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
  mono_indices <- which(maf == 0)
  if(length(mono_indices)==0){
    myGD_train = myGD_train
    myGD_test = myGD_test
  }else{
    myGD_train = myGD_train[,-mono_indices]
    myGD_test = myGD_test[,-mono_indices]}
  mydata <- data.frame(myY_train, myGD_train)
  colnames(mydata)[1]<-c("Y")
  gc()
  svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                   data = mydata,
                   method = model,
                   preProcess = c("center", "scale"),
                   trControl = trainControl(method = "cv", number = 10),
                   tuneLength = 10)
  svm.linear_pred1 <- predict(svmFit1, myGD_test)
  acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")

  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  SVM_acc_table <- rbind(SVM_acc_table, results_long)

  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}


NN_Mean_VS_M <- function(train_genotypes, train_phenotype,test_genotypes, test_phenotype,model="avNNet", folds = 5,markers=3000){
  gc()
  library(caret)
  library(tidyr)
  library(caret)    # for classification and regression training
  library(kernlab)
  library(randomForest)
  library(ranger)
  library(h2o)
  library(e1071)  # for fitting SVMs
  # Make the CV list
  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  myY_train <- train_phenotype
  myY_test <- test_phenotype
  samp=sample(1:length(train_genotypes), markers)

  myGD_train <- train_genotypes[,samp]
  myGD_test <- test_genotypes[,samp]

  maf <- calc_maf_apply(myGD_train, encoding = c(0, 1, 2))
  mono_indices <- which(maf == 0)
  #taxa <- row.names(myGD)
  if(length(mono_indices)==0){
    myGD_train = myGD_train
    myGD_test = myGD_test
  }else{
    myGD_train = myGD_train[,-mono_indices]
    myGD_test = myGD_test[,-mono_indices]}

  mydata <- data.frame(myY_train, myGD_train)
  colnames(mydata)[1]<-c("Y")
  gc()
  ## Create a specific candidate set of models to evaluate:
  nnetGrid <- expand.grid(decay = c(0, 0.01, .1),size = c(1:10),
                          ## The next option is to use bagging (see the
                          ## next chapter) instead of different random
                          ## seeds.
                          bag = FALSE)
  svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                   data = mydata,
                   method = model,
                   tuneGrid = nnetGrid,
                   trControl = trainControl(method = "cv", number = 10),
                   ## Automatically standardize data prior to modeling
                   ## and prediction
                   preProc = c("center", "scale"),
                   linout = TRUE,
                   trace = FALSE,
                   MaxNWts = 10 * (ncol(myGD_train) + 1) + 10 + 1,
                   maxit = 500)
  svm.linear_pred1 <- predict(svmFit1, myGD_test)
  acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")

  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  SVM_acc_table <- rbind(SVM_acc_table, results_long)

  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

NN_Mean_VS_PC <- function(train_PCA, train_phenotype,test_PCA,test_phenotype,model="avNNet", folds = 5)
{
  gc()
  library(caret)
  library(tidyr)
  library(caret)
  library(randomForest)# for classification and regression training
  library(kernlab)
  library(ranger)
  library(h2o)
  library(e1071)# for fitting SVMs
  # Make the CV list
  myY_train <- train_phenotype
  myY_test <- test_phenotype

  myPCA_train <- train_PCA[,c(1:ncol(test_PCA))]
  myPCA_test <- test_PCA

  gc()
  tooHigh <- findCorrelation(cor(myPCA_train), cutoff = .90)
  myPCA_train <- myPCA_train[-tooHigh,]
  myPCA_test <- myPCA_test[-tooHigh,]
  mydata <- data.frame(myY_train, myPCA_train)
  colnames(mydata)[1]<-c("Y")
  ## Create a specific candidate set of models to evaluate:
  nnetGrid <- expand.grid(decay = c(0, 0.01, .1),
                          size = c(1:10),
                          ## The next option is to use bagging (see the
                          ## next chapter) instead of different random
                          ## seeds.
                          bag = FALSE)
  svmFit1 <- train(y = mydata$Y, x = mydata[, colnames(mydata) != 'Y'],
                   data = mydata,
                   method = model,
                   tuneGrid = nnetGrid,
                   trControl = trainControl(method = "cv", number = 10),
                   ## Automatically standardize data prior to modeling
                   ## and prediction
                   preProc = c("center", "scale"),
                   linout = TRUE,
                   trace = FALSE,
                   MaxNWts = 10 * (ncol(myPCA_train) + 1) + 10 + 1,
                   maxit = 500)
  svm.linear_pred1 <- predict(svmFit1, myPCA_test)
  acc=cor(svm.linear_pred1, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")

  model_vect <- c(model)
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))

  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  SVM_acc_table <- rbind(SVM_acc_table, results_long)

  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}



BLUP_pc_mean_GAGS <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])

    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]

    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)

    GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
    if(length(GWASSM)==0){
      acc=NA
      results <- list(acc)
      names(results) <- c("acc")
    }else{
    sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

    myCV_train <- myGD_train[,sm]
    myCV_test <- myGD_test[,sm]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]

    fix_train <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test  <- as.matrix(cbind(myCV_test,myPCA_test))
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
    fix_train=fix_train[,-c(findLinearCombos(fix_train)$remove)]
    fix_test=fix_test[,-c(findLinearCombos(fix_test)$remove)]}
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    }
    BGLR_acc_results[[i]] <- list(acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

BLUP_pc_mean_GAGS_10 <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])

    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]

    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)

    top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
    GWASSM=top10[1:10,]$SNP
    myCV_train <- myGD_train[,GWASSM]
    myCV_test <- myGD_test[,GWASSM]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]

    fix_train <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test  <- as.matrix(cbind(myCV_test,myPCA_test))
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      fix_train=fix_train[,-c(findLinearCombos(fix_train)$remove)]
      fix_test=fix_test[,-c(findLinearCombos(fix_test)$remove)]}
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")

    BGLR_acc_results[[i]] <- list(acc)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

MAS_CV_Mean <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  svm_results <- list()
  length(fold_list)
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[1]])

    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_test <- as.matrix(genotypes[-fold_indices,])

    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]

    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)
    GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
    if(length(GWASSM)==0){
      acc=NA
      results <- list(acc)
      names(results) <- c("acc")
      svm_results[[i]] <- results
    }else{

      sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

      MAS_train <- myGD_train[,sm]
      MAS_test <- data.frame(myGD_test[,sm])
      str(MAS_test)
      gc()

      GLM_data <- data.frame(myY_train, MAS_train)
      str(GLM_data)
      names(GLM_data)[1] <- "Y"
      #Linear model to calculate effects
      #You can run all signficant markers at once to see cumulative effect
      MAS_model <- lm(Y ~ ., data = GLM_data)
      MAS_pred <- predict(MAS_model, MAS_test)
      acc=cor(MAS_pred, myY_test, use = "pairwise.complete")
      results <- list(acc)
      names(results) <- c("acc")
      svm_results[[i]] <- results
    }


  }
  model_vect <- c("MAS")
  SVM_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(svm_results)){
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(svm_results[[i]]))
    SVM_acc_table <- rbind(SVM_acc_table, results_long)
  }
  names(SVM_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(SVM_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

test_all_models_BLUP_vs_pc_mean_recode <- function(train_genotypes, train_phenotype,train_PCA=NULL,test_genotypes, test_phenotype,test_PCA=NULL)
{
  library(BGLR)
  library(tidyr)
  library(rrBLUP)
  library(caret)
  library(dplyr)

  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA

  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = myPCA_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- myPCA_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  metrics=postResample(pred=predictions,obs=myY_test)
  results=c(ACC=acc,metrics)
  prediction=data.frame(test_phenotype[,c(1,2)],GEBV=predictions,RE=pred_effects,FE=fix_effects)
  return(Accuracy=results,Predictions=prediction)
}


test_all_models_BLUP_vs_pc_mean_recode_M <- function(train_genotypes, train_phenotype,train_PCA=NULL,train_CV=NULL,test_genotypes, test_phenotype,test_PCA=NULL,test_CV=NULL)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA

  train_CV=impute(as.factor(train_CV))
  myCV_train=as.numeric(train_CV)

  test_CV=impute(as.factor(test_CV))
  myCV_test=as.numeric(test_CV)

  fix_train <- as.matrix(cbind(myCV_train,myPCA_train))
  fix_test  <- as.matrix(cbind(myCV_test,myPCA_test))
  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = fix_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- fix_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")


  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("VS", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

BLUP_vs_pc_mean_GAGS <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)


  # Split into training and testing data
  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  GD_train <- train_GD
  Y_train <-train_Y[,c(1,2)]

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]

  GWASR<- GAPIT(Y = Y_train,
                GD = GD_train,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)

  GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
  if(length(GWASSM)==0){
    acc=NA
    results <- list(acc)
    names(results) <- c("acc")
  }else{
    sm=as.character(GWASR$GWAS[GWASSM,]$SNP)

    myCV_train <- myGD_train[,sm]
    myCV_test <- myGD_test[,sm]
    myPCA_train <- train_PCA
    myPCA_test <- test_PCA

    fix_train <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test  <- as.matrix(cbind(myCV_test,myPCA_test))
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      fix_train=fix_train[,-c(findLinearCombos(fix_train)$remove)]
      fix_test=fix_test[,-c(findLinearCombos(fix_test)$remove)]}
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)

    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")

    results <- list(acc)
    names(results) <- c("acc")
  }
  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("VS", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]
  return(results)
}

BLUP_vs_pc_mean_GAGS_10 <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)

  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)

  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)

  GD_train <- train_GD
  Y_train <-train_Y[,c(1,2)]

  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]

  GWASR<- GAPIT(Y = Y_train,
                GD = GD_train,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)

  top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
  GWASSM=top10[1:10,]$SNP
  myCV_train <- myGD_train[,GWASSM]
  myCV_test <- myGD_test[,GWASSM]
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA

  fix_train <- as.matrix(cbind(myCV_train,myPCA_train))
  fix_test  <- as.matrix(cbind(myCV_test,myPCA_test))
  p <- ncol(fix_train)
  XtX <- crossprod(fix_train, fix_train)
  rank.X <- qr(XtX)$rank
  if (rank.X < p) {
    fix_train=fix_train[,-c(findLinearCombos(fix_train)$remove)]
    fix_test=fix_test[,-c(findLinearCombos(fix_test)$remove)]}
  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = fix_train)

  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- fix_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")

  results <- list(acc)
  names(results) <- c("acc")

  model_vect <- c("rrBLUP")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  results_long <- data.frame(rep(1, length(model_vect)), model_vect, unlist(results))
  BGLR_acc_table <- rbind(BGLR_acc_table, results_long)

  names(BGLR_acc_table) <- c("VS", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold)[2]

  return(results)
}

BGLR_cv_Ordinal_Mean <- function(genotypes, phenotype, nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA

    # Calculate the GS model using BGLR
    ##Ordinal
    BO_ETA<-list(list(X=as.matrix(genotypes),model="BL"))
    BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    BO_predictions <- predict(BO_results)
    BO_acc <- cor(as.numeric(phenotype[-fold_indices]), BO_predictions[-fold_indices])
    DF=BO_results$probs[-fold_indices,]
    probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
    probs=as.numeric(probs)
    BO_acc_cat=cor(as.numeric(phenotype[-fold_indices]), probs,use = "complete.obs")
    BGLR_acc_results[[i]] <- list(Acc=BO_acc,Acc_cat=BO_acc_cat)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("R_acc","Cat_acc")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:3])
  return(results)
}


BGLR_cv_Ordinal_pc_Mean <- function(genotypes, phenotype,PCA=NULL, nIter = 5000, burnIn = 2000, folds = 5)
{
  library(BGLR)
  library(tidyr)
  # Make the CV list


  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()

  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[pheno_train=="NaN"]<-NA
    pheno_train=droplevels(pheno_train)
    pheno_train[-fold_indices] <- NA

    # Calculate the GS model using BGLR
    ##Ordinal
    BO_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
    BO_results <- BGLR(y = pheno_train,response_type = 'ordinal', ETA = BO_ETA, nIter=nIter, burnIn=burnIn)
    BO_predictions <- predict(BO_results)
    BO_acc <- cor(as.numeric(phenotype[-fold_indices]), BO_predictions[-fold_indices],use = "complete.obs")
    DF=BO_results$probs[-fold_indices,]
    probs=colnames(DF)[max.col(replace(DF, cbind(seq_len(nrow(DF)), max.col(DF,ties.method="first")), -Inf), "first")]
    probs=as.numeric(probs)
    BO_acc_cat=cor(as.numeric(phenotype[-fold_indices]), probs,use = "complete.obs")
    BGLR_acc_results[[i]] <- list(Acc=BO_acc,Acc_cat=BO_acc_cat)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("R_acc","Cat_acc")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[,2:3])
  return(results)
}

GS_models_mean_pc_recode <- function(genotypes, phenotype,PCA=NULL,nIter = 5000, burnIn = 2000, folds = 5,model="rrBLUP")
{
  library(BGLR)
  library(tidyr)
  library(caret)
  library(caret)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype), k = folds)

  BGLR_acc_results <- list()
  Predictions_ALL<-c()
  model_list<-list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])

    if(model=="rrBLUP"){
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)

    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)

    PCA_bl_train=prcomp(myGD_train)
    myPCA_train=PCA_bl_train$x

    PCA_bl_test=prcomp(myGD_test)
    myPCA_test=PCA_bl_test$x

    myY_train <- phenotype[fold_indices]
    myY_test <- phenotype[-fold_indices]

    gc()
    model_results <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = myPCA_train)

    pred_effects <- myGD_test %*% model_results$u
    fix_effects <- myPCA_test  %*% model_results$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=myY_test,GEBV=predictions)
    }

    if(model=="GBLUP"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA

    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)

    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)

    PCA_bl_train=prcomp(myGD_train)
    myPCA_train=PCA_bl_train$x

    PCA_bl_test=prcomp(myGD_test)
    myPCA_test=PCA_bl_test$x
    gc()
    K <- tcrossprod(as.matrix(genotypes))
    #GBlup with fixed effects
    model_results <- mixed.solve(y = pheno_train,
                               K = K,
                               X=PCA)
    pred_effects <- model_results$u[-fold_indices]
    fix_effects <- as.matrix(myPCA_test)  %*% model_results$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=myY_test,GEBV=predictions)
    }

    if(model=="BayesA"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using BGLR
    ##Bayes A
    gc()
    bayesA_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesA"))
    model_results <- BGLR(y = pheno_train, ETA = bayesA_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('ba_',nIter,'_',burnIn,'_'))
    bayesA_predictions <- predict(model_results)
    acc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices], use = "pairwise.complete")
    sacc <- cor(phenotype[-fold_indices], bayesA_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=bayesA_predictions[-fold_indices],obs=phenotype[-fold_indices])
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=phenotype[-fold_indices],GEBV=bayesA_predictions[-fold_indices])
    }

    if(model=="BayesB"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    ##BayesB
    gc()
    bayesB_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesB"))
    model_results <- BGLR(y = pheno_train, ETA = bayesB_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('bb_',nIter,'_',burnIn,'_'))
    bayesB_predictions <- predict(model_results)
    acc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices], use = "pairwise.complete")
    sacc <- cor(phenotype[-fold_indices], bayesB_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=bayesB_predictions[-fold_indices],obs=phenotype[-fold_indices])
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=phenotype[-fold_indices],GEBV=bayesB_predictions[-fold_indices])
    }

    if(model=="BayesC"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    ##BayesC
    gc()
    bayesC_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BayesC"))
    model_results <- BGLR(y = pheno_train, ETA = bayesC_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('bc_',nIter,'_',burnIn,'_'))
    bayesC_predictions <- predict(model_results)
    acc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices], use = "pairwise.complete")
    sacc <- cor(phenotype[-fold_indices], bayesC_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=bayesC_predictions[-fold_indices],obs=phenotype[-fold_indices])
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=phenotype[-fold_indices],GEBV=bayesC_predictions[-fold_indices])
    }

    if(model=="BayesRR"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    ##BayesRR
    gc()
    BRR_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BRR"))
    model_results <- BGLR(y = pheno_train, ETA = BRR_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('brr_',nIter,'_',burnIn,'_'))
    BRR_predictions <- predict(model_results)
    acc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices], use = "pairwise.complete")
    sacc <- cor(phenotype[-fold_indices], BRR_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=BRR_predictions[-fold_indices],obs=phenotype[-fold_indices])
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=phenotype[-fold_indices],GEBV=BRR_predictions[-fold_indices])
    }

    if(model=="BayesL"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    ##Bayes L
    gc()
    BL_ETA<-list(list(X=PCA,model="FIXED"),list(X=as.matrix(genotypes),model="BL"))
    model_results <- BGLR(y = pheno_train, ETA = BL_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('bl_',nIter,'_',burnIn,'_'))
    BL_predictions <- predict(model_results)
    acc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices], use = "pairwise.complete")
    sacc <- cor(phenotype[-fold_indices], BL_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=BL_predictions[-fold_indices],obs=phenotype[-fold_indices])
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=phenotype[-fold_indices],GEBV=BL_predictions[-fold_indices])
    }

    if(model=="RKHS"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    ##Single Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    p<-ncol(X)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5
    K<-exp(-h*D)
    RKHS_ETA<-list(list(X=PCA,model="FIXED"),list(K=K,model='RKHS'))
    model_results <- BGLR(y = pheno_train, ETA = RKHS_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('RKHS_',nIter,'_',burnIn,'_'))
    RKHS_predictions <- predict(model_results)
    acc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices], use = "pairwise.complete")
    sacc <- cor(phenotype[-fold_indices], RKHS_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=RKHS_predictions[-fold_indices],obs=phenotype[-fold_indices])
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=phenotype[-fold_indices],GEBV=RKHS_predictions[-fold_indices])
    }

    if(model=="MK"){
    # Split into training and testing data
    pheno_train <- phenotype
    pheno_train[-fold_indices] <- NA
    ##Multi-Kernel RKHS
    gc()
    X=as.matrix(genotypes)
    X<-scale(X,center=TRUE,scale=TRUE)
    D<-(as.matrix(dist(X,method='euclidean'))^2)/p
    h<-0.5*c(1/5,1,5)
    MK_ETA<-list(list(X=PCA,model="FIXED"),
                 list(K=exp(-h[1]*D),model='RKHS'),
                 list(K=exp(-h[2]*D),model='RKHS'),
                 list(K=exp(-h[3]*D),model='RKHS'))
    model_results <- BGLR(y = pheno_train, ETA = MK_ETA, nIter=nIter, burnIn=burnIn,saveAt=paste0('MK_',nIter,'_',burnIn,'_'))
    MK_predictions <- predict(model_results)
    acc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices], use = "pairwise.complete")
    sacc <- cor(phenotype[-fold_indices], MK_predictions[-fold_indices], use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=MK_predictions[-fold_indices],obs=phenotype[-fold_indices])
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y=phenotype[-fold_indices],GEBV=MK_predictions[-fold_indices])
    }
    Predictions<-prediction
    Predictions_ALL=rbind(Predictions_ALL,Predictions)
    BGLR_acc_results[[i]] <- list(results)
    model_list[[i]]<-model_results


  }
  #, GBLUP_acc, "GBLUP"

  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold, na.rm = TRUE)[2:6]
  results_ALL=list(results,Predictions_ALL,model_list)
  return(results_ALL)
}

GE_Matrix <- function(genotypes, phenotype,trait,GE=FALSE){
  library(BGLR)
  library(BMTME)
  library(rrBLUP)
  library(tidyr)
  Pheno=phenotype
  Pheno=droplevels(Pheno)
  Pheno=complete(Pheno, Genotype, ENV)
if(GE==TRUE){
  LG <- cholesky(A.mat(genotypes))
  ZG <- model.matrix(~0 + as.factor(Pheno$Genotype))
  Z.G <- ZG %*% LG
  Z.E <- model.matrix(~0 + as.factor(Pheno$ENV))
  ZEG <- model.matrix(~0 + as.factor(Pheno$Genotype):as.factor(Pheno$ENV))
  G2 <- kronecker(diag(length(unique(Pheno$ENV))), data.matrix(A.mat(genotypes)))
  LG2 <- cholesky(G2)
  Z.EG <- ZEG %*% LG2
  Y <- as.matrix(Pheno[,trait])
  results=list(Y=Y,X=Z.E,Z1=Z.G,Z2=Z.EG)
}else{
  LG <- cholesky(A.mat(genotypes))
  ZG <- model.matrix(~0 + as.factor(Pheno$Genotype))
  Z.G <- ZG %*% LG
  Y <- as.matrix(Pheno[,trait])
  results=list(Y=Y,Z1=Z.G)
}

  return(results)
}
