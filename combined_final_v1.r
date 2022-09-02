library(caret)
library(corrplot)
library(dplyr)
mrna_dat <- read.table("C://users//aaswi//Downloads//MSKCC_PCa_mRNA_data (1).txt", header = TRUE)
mrna_updated=read.csv("C://users//aaswi//Downloads//mrna_updated.csv", header = TRUE)


mrna_t <- t(mrna_dat)

colnames(mrna_t)=mrna_t[2,]
mrna_t=mrna_t[-c(1,2),]

memory.limit(size = 56000)



mrna_t_updated=mrna_t

mrna_t.bkp=mrna_t
mrna_t_updated[,1:18564] <- lapply(mrna_t_updated[,1:18564], function(x) {
  if(is.character(x)) as.numeric(x) 
})




prostate_data  <- read.csv("C://users//aaswi//Downloads//data.csv")
prostate_data = prostate_data[,c(1,23)]
prostate_data=na.omit(prostate_data)

prostate_data$PathStage<-
  plyr::mapvalues(prostate_data$PathStage, from = c("T2A","T2B", "T2C","T3A", "T3B", "T3C","T4"),
                  to =c("Early", "Early", "Early","Advanced", "Advanced", "Advanced","Advanced"))

path_data <- prostate_data
rownames(path_data)=path_data[,1]
#path_data=path_data[,-c(1)]



comb_data=merge(path_data,mrna_t_updated,by=0)
rownames(comb_data)=comb_data[,1]
comb_data=comb_data[,-c(1,2)]
#comb_data=comb_data[,-c("PathStage")]
comb_data=subset(comb_data,select = -c(PathStage))

##############clinical
clin = read.csv("C://users//aaswi//Downloads//data.csv")
dat = clin

dat=dat[-c(27:39)]
colnames(dat)

#dat.merge=merge(dat,var_df,by="Sample.ID")
#dat=dat.merge
#dat.t=t(dat)


rownames(dat)=dat[,1]
dat=dat[-c(1)]

library(plyr)
library(dplyr)
## Remove columns with more than 50% NA

dat=dat[which(rowMeans(!is.na(dat)) > 0.5), ]

dat <- subset(dat, !is.na(dat$PathStage))
dat <- subset(dat, !is.na(dat$Race))
is.na(dat$ClinT_Stage)
dat <- subset(dat, !is.na(dat$ClinT_Stage))


dat$Race <-
  mapvalues(dat$Race, from = c("Black Non Hispanic", "White Non Hispanic",
                               "Black Hispanic", "White Hispanic", "Asian","Unknown"),
            to =c("Black", "White", "Black", "White", "Asian","Unknown"))

dat$ClinT_Stage <-
  mapvalues(dat$ClinT_Stage, from = c("T1C", "T2", "T2A","T2B", "T2C", "T3",
                                      "T3A", "T3B", "T3C"),
            to =c("T1", "T2", "T2","T2", "T2", "T3","T3", "T3", "T3"))
dat$PathStage<-
  mapvalues(dat$PathStage, from = c("T2A","T2B", "T2C","T3A", "T3B", "T3C","T4"),
            to =c("Early", "Early", "Early","Advanced", "Advanced", "Advanced","Advanced"))


#unique(dat$MetSite)

#?subset()

colnames(dat)

#nrow(dat)
#is.na(dat$MetSite)

dat$MetSite <-
  mapvalues(dat$MetSite, from = c(NA, "spine",
                                  "bone", "brain", "exent bladder locally progressive","bone lung","node","exent colon locally progressive"),
            to =c("No", "spine",
                  "bone", "brain", "exent bladder locally progressive","bone lung","node","exent colon locally progressive"))

dat <- subset(dat, !is.na(dat$MetSite))

#is.na(dat$PreDxBxPSA)
dat <- subset(dat, !is.na(dat$PreDxBxPSA))




#is.na(dat$PathGG1)

dat$NeoAdjRadTx <-
  mapvalues(dat$NeoAdjRadTx, from = c(NA),
            to =c("No"))

dat$ChemoTx <-
  mapvalues(dat$ChemoTx, from = c(NA),
            to =c("No"))

dat$HormTx <-
  mapvalues(dat$HormTx, from = c(NA),
            to =c("No"))

dat$RadTxType <-
  mapvalues(dat$RadTxType, from = c(NA),
            to =c("No"))

####remove num nodes removed and positive

colnames(dat)
dat=dat[-c(20:21)]



dat <- subset(dat, !is.na(dat$PathGG1))
dat <- subset(dat, !is.na(dat$PathGG2))
dat <- subset(dat, !is.na(dat$PathGGS))


#unique(dat$PathGGS)

colnames(dat)


library(dplyr)
dat <- dat %>% 
  mutate(Tx = if_else(NeoAdjRadTx != "No" | ChemoTx != "No" |
                        HormTx != "No" |RadTxType !="No", "Therapy_done", "No therapy"))
colnames(dat)
#dat[c(11,12,13,14,24)]

#nrow(dat)




#remove Tx columns
dat=dat[-c(11:14)]



#dat["Type"]



#dat[rowSums(is.na(dat)) > 0, ]     


#class(dat$Pathstage)





##########rfe and randomforest
dat$Type=as.factor(dat$Type)
dat$MetSite=as.factor(dat$MetSite)
dat$Race=as.factor(dat$Race)
dat$ClinT_Stage=as.factor(dat$ClinT_Stage)
dat$RP_Type=as.factor(dat$RP_Type)
dat$SMS=as.factor(dat$SMS)
dat$ECE=as.factor(dat$ECE)
dat$SVI=as.factor(dat$SVI)
dat$LNI=as.factor(dat$LNI)
dat$PathStage=as.factor(dat$PathStage)
dat$BxGG1=as.factor(dat$BxGG1)
dat$BxGG2=as.factor(dat$BxGG2)
dat$BxGGS=as.factor(dat$BxGGS)
dat$PathGG1=as.factor(dat$PathGG1)
dat$PathGG2=as.factor(dat$PathGG2)
dat$PathGGS=as.factor(dat$PathGGS)
dat$Tx=as.factor(dat$Tx)

dat1=merge(dat,comb_data,by=0)
#dat1=merge(dat,mrna_t_updated,by=0)
colnames(dat1)=str_replace_all(colnames(dat1), "[^[:alnum:]]", "")
dat1=dat1[,-c(1,2,5,12,13,14,15,16,18,19,20)]
pthstg=dat1$PathStage
dat1=dat1[,-9]
dat1$PathStage=pthstg
dat1$PathStage=as.factor(dat1$PathStage)
#names(dat1) <- gsub(" ", ".", names(dat1))
pros_new_dat = dat1
R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
sensitivities = c()
specificities = c()
thresholds = c()

#par(mfrow=c(3,3))
set.seed(6090)
library(pROC)
library(caret)
library(randomForest)
library(caret)
library(leaps)
#vec=seq(1,55)
for(r in 1:5){
  n = nrow(pros_new_dat)
  i.cv = sample(1:n,n, replace=FALSE)
  print(i.cv[1])
  pc_copy = pros_new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.rf = numeric(K)
  feature.matrix = matrix(0,nrow = K, ncol = ncol(pc_train))
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-18574]
  
  #consolidated_list = list()
  
  
  #set.seed(6090)
  for(k in 1:K){
    
    i.train	= which(folds!=k)
    pros_train_old = pc_train[i.train,]
    pros_valid_old = pc_train[-i.train,]
    pc_numeric_pros_train = pros_train_old[,-c(1,2,4,5,6,8,9)]
    pc_categorical_pros_train = pros_train_old[,c(1,2,4,5,6,8,9,18574)]
    
    #chi sqr
    features=c()
    categorical_features=c()
    temp=c()
    adjust=c()
    for(col in 1:(ncol(pc_categorical_pros_train)-1)){
      chisq_test <- chisq.test(pc_categorical_pros_train[,col], pc_categorical_pros_train$PathStage, correct = FALSE)
      temp[col]=chisq_test$p.value
      #if(chisq_test$p.value<0.05){
      #  categorical_features = c(categorical_features,colnames(pc_categorical_pros_train)[col])
      #}
    }
    adjust=p.adjust(temp,method = "fdr")
    for(j in 1:(ncol(pc_categorical_pros_train)-1)){
      if(adjust[j]<0.05){
        categorical_features=c(categorical_features,colnames(pc_categorical_pros_train)[j])
      }
    }
    #mutual info
    numerical_features=c()
    mrna.feature.matrix.mmi = matrix(nrow = 1, ncol = ncol(pc_numeric_pros_train)-1)
    colnames(mrna.feature.matrix.mmi) = colnames(pc_numeric_pros_train[1:18566])
    for(i in 1:18566){
      mmi.val=mmi(as.matrix(pc_numeric_pros_train[,i]),data.frame(pc_numeric_pros_train$PathStage))
      mrna.feature.matrix.mmi[1,i]=mmi.val$mi
    }
    mmi.val09.df=data.frame(mrna.feature.matrix.mmi)
    mmi.val09.df=mmi.val09.df[colSums(mmi.val09.df>0.13)>0]
    numerical_features=c(numerical_features,colnames(mmi.val09.df))
    
    
    features=c(numerical_features,categorical_features,"PathStage")
    cat("r = ",r,"k = ",k,"sel features = ",features)
    
    pros_train=pros_train_old[(colnames(pros_train_old) %in% features)]
    pros_valid=pros_valid_old[(colnames(pros_valid_old) %in% features)]
    
    #pros_train=pros_train_old[,c(features)]
    #pros_valid=pros_valid_old[,c(features)]
    
    x.train = pros_train[,-c(length(features))]
    y.train = pros_train$PathStage
    x.valid = pros_valid[,-c(length(features))]
    y.valid = pros_valid$PathStage
    #subsets <- c(1,2,3, ncol(pros_train[,-c(11)]))
    vec=seq(1,length(features))
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "cv",
                       number = 10,
                       # method = "repeatedcv",
                       # repeats = 5,
                       verbose = FALSE)
    rf.rfe <- rfe(x.train, y.train,
                  sizes = vec,
                  #sizes=subsets,
                  rfeControl = ctrl)
    vars = match(rf.rfe$optVariables,colnames(pros_train))
    feature.matrix[k,rf.rfe$optVariables] = 1
    #feature.matrix[k,match(rf.rfe$optVariables,colnames(x.train))] = 1
    #feature.matrix[k,-match(rf.rfe$optVariables,colnames(x.train))] = 0
    pros_fs_train = pros_train[,c(vars,match('PathStage',colnames(pros_train)))]
    pros_fs_valid = pros_valid[,c(vars,match('PathStage',colnames(pros_valid)))]
    
    
    #pros_fs_train = pros_train[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_train)))]
    #pros_fs_valid = pros_valid[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_valid)))]
    rf.o = randomForest(PathStage~., data=pros_fs_train)
    rf.p = predict(rf.o, newdata=pros_fs_valid)
    tb.rf = table(rf.p, pros_fs_valid$PathStage)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    #consolidated_list <- append(consolidated_list, list(c(vars,acc.rf[k])))
    cat("r = ",r,"k = ",k,"internal acc = ",acc.rf[k])
  }
  #print(consolidated_list)
  accuracy.avg.inner[r] = mean(acc.rf)
  cat("r = ",r,"mean accuracy is :",accuracy.avg.inner[r])
  cat("r = ",r,"SD of accuracy is :",sd(acc.rf))
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.5){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  print(fs)
  #fs=unique(fs)
  #print(fs)
  pc_train_fs = pc_train[,fs] 
  pc_test_fs = pc_test[,fs] 
  rf.fit = randomForest(PathStage~., data=pc_train_fs)
  rf.pred = predict(rf.fit, newdata=pc_test_fs)
  table.rf = confusionMatrix(rf.pred, pc_test_fs$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.rf$overall[1],2)
  precision[r] = round(table.rf$byClass[5],2)
  recall[r] = round(table.rf$byClass[6],2)
  f1.score[r] = round(table.rf$byClass[7],2)
  rf.p = predict(rf.fit, pc_test_fs, type="prob")[,2] 
  roc.rf = roc(response=pc_test_fs$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.rf
  auc.outer[r] = round(roc.rf$auc,2)
  sensitivities[r] = list(roc.rf$sensitivities)
  specificities[r] = list(roc.rf$specificities)
  thresholds[r] = list(roc.rf$thresholds)
}

plot(roc.outer[[1]], col='black')
plot(roc.outer[[2]], add=TRUE, col='red')
plot(roc.outer[[3]], add=TRUE, col='skyblue')
plot(roc.outer[[4]], add=TRUE, col='green')
plot(roc.outer[[5]], add=TRUE, col='orange')
legend(x = 'bottomright', legend=c("Iter 1", "Iter 2", "Iter 3","Iter 4","Iter 5"), 
       col = c("black","red",'skyblue', 'green','orange'), lty=1, cex=1.0, title = 'RandomForest'
)

max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  #mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i])/2
  
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  #mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i])/2
  
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities
par(pty='s')
my_labels <- sprintf(rev(seq(0.0, 1.0, 0.2)),
                     fmt = '%#.1f')


par(pty='s')
plot(1-new_specificities, new_sensitivities, type='l',xlab='1-Specificity',ylab='Sensitivity')
abline(coef=c(0,1))

#############stepwise and SVM
dat$Type=as.factor(dat$Type)
dat$MetSite=as.factor(dat$MetSite)
dat$Race=as.factor(dat$Race)
dat$ClinT_Stage=as.factor(dat$ClinT_Stage)
dat$RP_Type=as.factor(dat$RP_Type)
dat$SMS=as.factor(dat$SMS)
dat$ECE=as.factor(dat$ECE)
dat$SVI=as.factor(dat$SVI)
dat$LNI=as.factor(dat$LNI)
dat$PathStage=as.factor(dat$PathStage)
dat$BxGG1=as.factor(dat$BxGG1)
dat$BxGG2=as.factor(dat$BxGG2)
dat$BxGGS=as.factor(dat$BxGGS)
dat$PathGG1=as.factor(dat$PathGG1)
dat$PathGG2=as.factor(dat$PathGG2)
dat$PathGGS=as.factor(dat$PathGGS)
dat$Tx=as.factor(dat$Tx)


dat$Type=as.numeric(dat$Type)
dat$MetSite=as.numeric(dat$MetSite)
dat$Race=as.numeric(dat$Race)
dat$ClinT_Stage=as.numeric(dat$ClinT_Stage)
dat$RP_Type=as.numeric(dat$RP_Type)
dat$SMS=as.numeric(dat$SMS)
dat$ECE=as.numeric(dat$ECE)
dat$SVI=as.numeric(dat$SVI)
dat$LNI=as.numeric(dat$LNI)
dat$BxGG1=as.numeric(dat$BxGG1)
dat$BxGG2=as.numeric(dat$BxGG2)
dat$BxGGS=as.numeric(dat$BxGGS)
dat$PathGG1=as.numeric(dat$PathGG1)
dat$PathGG2=as.numeric(dat$PathGG2)
dat$PathGGS=as.numeric(dat$PathGGS)
dat$Tx=as.numeric(dat$Tx)

dat1=merge(dat,comb_data,by=0)
colnames(dat1)=str_replace_all(colnames(dat1), "[^[:alnum:]]", "")
dat1=dat1[,-c(1,2,5,12,13,14,15,16,18,19,20)]
pthstg=dat1$PathStage
dat1=dat1[,-9]
dat1$PathStage=pthstg
dat1$PathStage=as.factor(dat1$PathStage)
#names(dat1) <- gsub(" ", ".", names(dat1))
pros_new_dat = dat1
R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
sensitivities = c()
specificities = c()
thresholds = c()

#par(mfrow=c(3,3))
set.seed(6090)
library(pROC)
library(caret)
library(randomForest)
library(caret)
library(leaps)
library(e1071)

for(r in 1:5){
  n = nrow(pros_new_dat)
  i.cv = sample(1:n, replace=FALSE)
  print(i.cv[1])
  pc_copy = pros_new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.svm = numeric(K)
  feature.matrix = matrix(0,nrow = K, ncol = ncol(pc_train))
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-18574]
  
  consolidated_list = list()
  
  #set.seed(6090)
  for(k in 1:10){
    i.train	= which(folds!=k)
    #pros_train = pc_train[i.train,]
    #pros_valid = pc_train[-i.train,]
    
    pros_train_old = pc_train[i.train,]
    pros_valid_old = pc_train[-i.train,]
    pc_numeric_pros_train = pros_train_old[,-c(1,2,4,5,6,8,9)]
    pc_categorical_pros_train = pros_train_old[,c(1,2,4,5,6,8,9,18574)]
    
    #chi sqr
    features=c()
    categorical_features=c()
    temp=c()
    adjust=c()
    for(col in 1:(ncol(pc_categorical_pros_train)-1)){
      chisq_test <- chisq.test(pc_categorical_pros_train[,col], pc_categorical_pros_train$PathStage, correct = FALSE)
      temp[col]=chisq_test$p.value
      #if(chisq_test$p.value<0.05){
      #  categorical_features = c(categorical_features,colnames(pc_categorical_pros_train)[col])
      #}
    }
    adjust=p.adjust(temp,method = "fdr")
    for(j in 1:(ncol(pc_categorical_pros_train)-1)){
      if(adjust[j]<0.05){
        categorical_features=c(categorical_features,colnames(pc_categorical_pros_train)[j])
      }
    }
    
    #mutual info
    numerical_features=c()
    mrna.feature.matrix.mmi = matrix(nrow = 1, ncol = ncol(pc_numeric_pros_train)-1)
    colnames(mrna.feature.matrix.mmi) = colnames(pc_numeric_pros_train[1:18566])
    for(i in 1:18566){
      mmi.val=mmi(as.matrix(pc_numeric_pros_train[,i]),data.frame(pc_numeric_pros_train$PathStage))
      mrna.feature.matrix.mmi[1,i]=mmi.val$mi
    }
    mmi.val09.df=data.frame(mrna.feature.matrix.mmi)
    mmi.val09.df=mmi.val09.df[colSums(mmi.val09.df>0.13)>0]
    numerical_features=c(numerical_features,colnames(mmi.val09.df))
    
    
    features=c(numerical_features,categorical_features,"PathStage")
    cat("r = ",r,"k = ",k,"sel features = ",features)
    
    pros_train=pros_train_old[(colnames(pros_train_old) %in% features)]
    pros_valid=pros_valid_old[(colnames(pros_valid_old) %in% features)]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    feature.matrix[k,names(which(mod)[-1]-1)] = 1
    #feature.matrix[k,-names(which(mod)[-1]-1)] = 0
    vars = rownames(data.frame(which(mod)[-1]-1))
    pros_fs_train = pros_train[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_train)))]
    pros_fs_valid = pros_valid[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_valid)))]
    
    cv.xtrain=pros_fs_train[, !colnames(pros_fs_train) %in% c("PathStage")]
    cv.ytrain=pros_fs_train[, colnames(pros_fs_train) %in% c("PathStage")]
    cv.xtest=pros_fs_valid[, !colnames(pros_fs_valid) %in% c("PathStage")]
    cv.ytest=pros_fs_valid[, colnames(pros_fs_valid) %in% c("PathStage")]
    
    xm = as.matrix(cv.xtrain)
    svm.tune = e1071::tune(svm, train.x=cv.xtrain, train.y=cv.ytrain,
                           kernel='polynomial',
                           ranges=list(cost=10^(-2:4), 
                                       #ranges=list(cost=10,
                                       #gamma=c(0.25)))
                                       gamma=c(0.25,0.5,1,1.5,2)))
    svmo.final = svm(xm, cv.ytrain, kernel='polynomial',
                     gamma=svm.tune$best.parameters$gamma,
                     cost=svm.tune$best.parameters$cost)
    pre=predict(svmo.final,cv.xtest)
    tb.svm=table(pre, cv.ytest)
    
    #rf.o = randomForest(PathStage~., data=pros_fs_train)
    #rf.p = predict(rf.o, newdata=pros_fs_valid)
    #tb.rf = table(rf.p, pros_fs_valid$PathStage)
    acc.svm[k] = sum(diag(tb.svm)) / sum(tb.svm)
    #consolidated_list <- append(consolidated_list, list(c(vars,acc.svm[k])))
    cat("r = ",r,"k = ",k,"internal acc = ",acc.svm[k])
  }
  #print(consolidated_list)
  accuracy.avg.inner[r] = mean(acc.svm)
  cat("r = ",r,"mean accuracy is :",accuracy.avg.inner[r])
  cat("r = ",r,"SD of accuracy is :",sd(acc.svm))
  
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.5){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  print(fs)
  #fs=unique(fs)
  #print(fs)
  pc_train_fs = pc_train[,fs] 
  pc_test_fs = pc_test[,fs] 
  
  outer.cv.xtrain=pc_train_fs[, !colnames(pc_train_fs) %in% c("PathStage")]
  outer.cv.ytrain=pc_train_fs[, colnames(pc_train_fs) %in% c("PathStage")]
  outer.cv.xtest=pc_test_fs[, !colnames(pc_test_fs) %in% c("PathStage")]
  outer.cv.ytest=pc_test_fs[, colnames(pc_test_fs) %in% c("PathStage")]
  
  xm = as.matrix(outer.cv.xtrain)
  svm.tune = e1071::tune(svm, train.x=outer.cv.xtrain, train.y=outer.cv.ytrain,
                         kernel='polynomial',
                         ranges=list(cost=10^(-2:4), 
                                     #ranges=list(cost=10,
                                     #gamma=c(0.25))
                                     gamma=c(0.25,0.5,1,1.5,2)))
  svmo.final = svm(xm, outer.cv.ytrain, kernel='polynomial',
                   gamma=svm.tune$best.parameters$gamma,
                   cost=svm.tune$best.parameters$cost)
  outer.pre=predict(svmo.final,outer.cv.xtest)
  
  
  #rf.fit = randomForest(PathStage~., data=pc_train_fs)
  #rf.pred = predict(rf.fit, newdata=pc_test_fs)
  #par(3,2)
  table.svm = confusionMatrix(outer.pre, outer.cv.ytest, mode = 'everything')
  accuracy.outer[r] = round(table.svm$overall[1],2)
  precision[r] = round(table.svm$byClass[5],2)
  recall[r] = round(table.svm$byClass[6],2)
  f1.score[r] = round(table.svm$byClass[7],2)
  svmo.roc= svm(xm, outer.cv.ytrain, kernel='polynomial',
                gamma=svm.tune$best.parameters$gamma,
                cost=svm.tune$best.parameters$cost,probability = TRUE)
  svm.p = predict(svmo.roc, outer.cv.xtest, probability = TRUE)
  attrb=attr(svm.p,"probabilities")[,2]
  roc.svm = roc(response=outer.cv.ytest, predictor=attrb)
  roc.outer[[r]] = roc.svm
  auc.outer[r] = round(roc.svm$auc,2)
  precision[r] = round(table.svm$byClass[5],2)
  recall[r] = round(table.svm$byClass[6],2)
  f1.score[r] = round(table.svm$byClass[7],2)
  sensitivities[r] = list(roc.svm$sensitivities)
  specificities[r] = list(roc.svm$specificities)
  thresholds[r] = list(roc.svm$thresholds)
  
}
plot(roc.outer[[1]], col='black')
plot(roc.outer[[2]], add=TRUE, col='red')
plot(roc.outer[[3]], add=TRUE, col='skyblue')
plot(roc.outer[[4]], add=TRUE, col='green')
plot(roc.outer[[5]], add=TRUE, col='orange')
legend(x = 'bottomright', legend=c("Iter 1", "Iter 2", "Iter 3","Iter 4","Iter 5"), 
       col = c("black","red",'skyblue', 'green','orange'), lty=1, cex=1.0, title = 'SVM'
)

max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  #mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i])/2
  
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  #mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i])/2
  
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities
par(pty='s')
my_labels <- sprintf(rev(seq(0.0, 1.0, 0.2)),
                     fmt = '%#.1f')


par(pty='s')
plot(1-new_specificities, new_sensitivities, type='l',xlab='1-Specificity',ylab='Sensitivity')
abline(coef=c(0,1))


##########stepwise and Decision tree
dat$Type=as.factor(dat$Type)
dat$MetSite=as.factor(dat$MetSite)
dat$Race=as.factor(dat$Race)
dat$ClinT_Stage=as.factor(dat$ClinT_Stage)
dat$RP_Type=as.factor(dat$RP_Type)
dat$SMS=as.factor(dat$SMS)
dat$ECE=as.factor(dat$ECE)
dat$SVI=as.factor(dat$SVI)
dat$LNI=as.factor(dat$LNI)
dat$PathStage=as.factor(dat$PathStage)
dat$BxGG1=as.factor(dat$BxGG1)
dat$BxGG2=as.factor(dat$BxGG2)
dat$BxGGS=as.factor(dat$BxGGS)
dat$PathGG1=as.factor(dat$PathGG1)
dat$PathGG2=as.factor(dat$PathGG2)
dat$PathGGS=as.factor(dat$PathGGS)
dat$Tx=as.factor(dat$Tx)


dat$Type=as.numeric(dat$Type)
dat$MetSite=as.numeric(dat$MetSite)
dat$Race=as.numeric(dat$Race)
dat$ClinT_Stage=as.numeric(dat$ClinT_Stage)
dat$RP_Type=as.numeric(dat$RP_Type)
dat$SMS=as.numeric(dat$SMS)
dat$ECE=as.numeric(dat$ECE)
dat$SVI=as.numeric(dat$SVI)
dat$LNI=as.numeric(dat$LNI)
dat$BxGG1=as.numeric(dat$BxGG1)
dat$BxGG2=as.numeric(dat$BxGG2)
dat$BxGGS=as.numeric(dat$BxGGS)
dat$PathGG1=as.numeric(dat$PathGG1)
dat$PathGG2=as.numeric(dat$PathGG2)
dat$PathGGS=as.numeric(dat$PathGGS)
dat$Tx=as.numeric(dat$Tx)

dat1=merge(dat,comb_data,by=0)
colnames(dat1)=str_replace_all(colnames(dat1), "[^[:alnum:]]", "")
dat1=dat1[,-c(1,2,5,12,13,14,15,16,18,19,20)]
pthstg=dat1$PathStage
dat1=dat1[,-9]
dat1$PathStage=pthstg
dat1$PathStage=as.factor(dat1$PathStage)
#names(dat1) <- gsub(" ", ".", names(dat1))
pros_new_dat = dat1
R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
sensitivities = c()
specificities = c()
thresholds = c()

#par(mfrow=c(3,3))
set.seed(6090)
library(pROC)
library(caret)
library(randomForest)
library(caret)
library(leaps)
library(e1071)

for(r in 1:5){
  n = nrow(pros_new_dat)
  i.cv = sample(1:n,n, replace=FALSE)
  print(i.cv[1])
  pc_copy = pros_new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  
  
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.dt = numeric(K)
  feature.matrix = matrix(0,nrow = K, ncol = ncol(pc_train))
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-18574]
  
  #consolidated_list = list()
  
  #set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    
    pros_train_old = pc_train[i.train,]
    pros_valid_old = pc_train[-i.train,]
    pc_numeric_pros_train = pros_train_old[,-c(1,2,4,5,6,8,9)]
    pc_categorical_pros_train = pros_train_old[,c(1,2,4,5,6,8,9,18574)]
    
    #chi sqr
    features=c()
    categorical_features=c()
    temp=c()
    adjust=c()
    for(col in 1:(ncol(pc_categorical_pros_train)-1)){
      chisq_test <- chisq.test(pc_categorical_pros_train[,col], pc_categorical_pros_train$PathStage, correct = FALSE)
      temp[col]=chisq_test$p.value
      #if(chisq_test$p.value<0.05){
      #  categorical_features = c(categorical_features,colnames(pc_categorical_pros_train)[col])
      #}
    }
    adjust=p.adjust(temp,method = "fdr")
    for(j in 1:(ncol(pc_categorical_pros_train)-1)){
      if(adjust[j]<0.05){
        categorical_features=c(categorical_features,colnames(pc_categorical_pros_train)[j])
      }
    }
    
    #mutual info
    numerical_features=c()
    mrna.feature.matrix.mmi = matrix(nrow = 1, ncol = ncol(pc_numeric_pros_train)-1)
    colnames(mrna.feature.matrix.mmi) = colnames(pc_numeric_pros_train[1:18566])
    for(i in 1:18566){
      mmi.val=mmi(as.matrix(pc_numeric_pros_train[,i]),data.frame(pc_numeric_pros_train$PathStage))
      mrna.feature.matrix.mmi[1,i]=mmi.val$mi
    }
    mmi.val09.df=data.frame(mrna.feature.matrix.mmi)
    mmi.val09.df=mmi.val09.df[colSums(mmi.val09.df>0.13)>0]
    numerical_features=c(numerical_features,colnames(mmi.val09.df))
    
    
    features=c(numerical_features,categorical_features,"PathStage")
    cat("r = ",r,"k = ",k,"sel features = ",features)
    
    pros_train=pros_train_old[(colnames(pros_train_old) %in% features)]
    pros_valid=pros_valid_old[(colnames(pros_valid_old) %in% features)]
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    feature.matrix[k,names(which(mod)[-1]-1)] = 1
    #feature.matrix[k,-(which(mod)[-1]-1)] = 0
    vars = rownames(data.frame(which(mod)[-1]-1))
    pros_fs_train = pros_train[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_train)))]
    pros_fs_valid = pros_valid[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_valid)))]
    dt.o = tree(PathStage~., data=pros_fs_train)
    dt.p = predict(dt.o, newdata=pros_fs_valid,type="class")
    tb.dt = table(dt.p, pros_fs_valid$PathStage)
    acc.dt[k] = sum(diag(tb.dt)) / sum(tb.dt)
    #consolidated_list <- append(consolidated_list, list(c(vars,acc.dt[k])))
    cat("r = ",r,"k = ",k,"internal acc = ",acc.dt[k])
  }
  #print(consolidated_list)
  accuracy.avg.inner[r] = mean(acc.dt)
  cat("r = ",r,"mean accuracy is :",accuracy.avg.inner[r])
  cat("r = ",r,"SD of accuracy is :",sd(acc.dt))
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.5){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  print(fs)
  #fs=unique(fs)
  #print(fs)
  pc_train_fs = pc_train[,fs] 
  pc_test_fs = pc_test[,fs] 
  dt.fit = tree(PathStage~., data=pc_train_fs)
  dt.pred = predict(dt.fit, newdata=pc_test_fs,type="class")
  #acc.dt = sum(diag(tb.dt)) / sum(tb.dt)
  table.dt = confusionMatrix(dt.pred, pc_test_fs$PathStage, mode = 'everything')
  accuracy.outer[r] = round(table.dt$overall[1],2)
  
  cv.CS = cv.tree(dt.fit, FUN=prune.misclass)
  opt.size = cv.CS$size[which.min(cv.CS$dev)]
  ptree = prune.misclass(dt.fit, best=opt.size)
  ptree.summary=summary(ptree)
  miss.class.err = 1-(ptree.summary$misclass[1]/ptree.summary$misclass[2])
  cat("r = ",r,"correct classification rate =",miss.class.err)
  print(ptree.summary$used)
  #roc.pred=as.numeric(predict(ptree,pc_test_fs,type='class'))
  
  precision[r] = round(table.dt$byClass[5],2)
  recall[r] = round(table.dt$byClass[6],2)
  f1.score[r] = round(table.dt$byClass[7],2)
  dt.p = predict(dt.fit, pc_test_fs, type="vector",probability=TRUE)[,2]
  roc.dt = roc(response=pc_test_fs$PathStage, predictor=dt.p)
  roc.outer[[r]] = roc.dt
  auc.outer[r] = round(roc.dt$auc,2)
  sensitivities[r] = list(roc.dt$sensitivities)
  specificities[r] = list(roc.dt$specificities)
  thresholds[r] = list(roc.dt$thresholds)
}

plot(roc.outer[[1]], col='black')
plot(roc.outer[[2]], add=TRUE, col='red')
plot(roc.outer[[3]], add=TRUE, col='skyblue')
plot(roc.outer[[4]], add=TRUE, col='green')
plot(roc.outer[[5]], add=TRUE, col='orange')
legend(x = 'bottomright', legend=c("Iter 1", "Iter 2", "Iter 3","Iter 4","Iter 5"), 
       col = c("black","red",'skyblue', 'green','orange'), lty=1, cex=1.0, title = 'Tree'
)

max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  #mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i])/2
  
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  #mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i])/2
  
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities
par(pty='s')
my_labels <- sprintf(rev(seq(0.0, 1.0, 0.2)),
                     fmt = '%#.1f')


par(pty='s')
plot(1-new_specificities, new_sensitivities, type='l',xlab='1-Specificity',ylab='Sensitivity')
abline(coef=c(0,1))


##########stepwise and Logistic Regression
dat$Type=as.factor(dat$Type)
dat$MetSite=as.factor(dat$MetSite)
dat$Race=as.factor(dat$Race)
dat$ClinT_Stage=as.factor(dat$ClinT_Stage)
dat$RP_Type=as.factor(dat$RP_Type)
dat$SMS=as.factor(dat$SMS)
dat$ECE=as.factor(dat$ECE)
dat$SVI=as.factor(dat$SVI)
dat$LNI=as.factor(dat$LNI)
dat$PathStage=as.factor(dat$PathStage)
dat$BxGG1=as.factor(dat$BxGG1)
dat$BxGG2=as.factor(dat$BxGG2)
dat$BxGGS=as.factor(dat$BxGGS)
dat$PathGG1=as.factor(dat$PathGG1)
dat$PathGG2=as.factor(dat$PathGG2)
dat$PathGGS=as.factor(dat$PathGGS)
dat$Tx=as.factor(dat$Tx)


dat$Type=as.numeric(dat$Type)
dat$MetSite=as.numeric(dat$MetSite)
dat$Race=as.numeric(dat$Race)
dat$ClinT_Stage=as.numeric(dat$ClinT_Stage)
dat$RP_Type=as.numeric(dat$RP_Type)
dat$SMS=as.numeric(dat$SMS)
dat$ECE=as.numeric(dat$ECE)
dat$SVI=as.numeric(dat$SVI)
dat$LNI=as.numeric(dat$LNI)
dat$BxGG1=as.numeric(dat$BxGG1)
dat$BxGG2=as.numeric(dat$BxGG2)
dat$BxGGS=as.numeric(dat$BxGGS)
dat$PathGG1=as.numeric(dat$PathGG1)
dat$PathGG2=as.numeric(dat$PathGG2)
dat$PathGGS=as.numeric(dat$PathGGS)
dat$Tx=as.numeric(dat$Tx)

dat1=merge(dat,comb_data,by=0)
colnames(dat1)=str_replace_all(colnames(dat1), "[^[:alnum:]]", "")
dat1=dat1[,-c(1,2,5,12,13,14,15,16,18,19,20)]
pthstg=dat1$PathStage
dat1=dat1[,-9]
dat1$PathStage=pthstg
dat1$PathStage=as.factor(dat1$PathStage)
#names(dat1) <- gsub(" ", ".", names(dat1))
pros_new_dat = dat1
R = 5
accuracy.outer = numeric(R)
accuracy.avg.inner = numeric(R)
auc.outer = numeric(R)
roc.outer = list(R)
precision = numeric(R)
recall = numeric(R)
f1.score = numeric(R)
sensitivities = c()
specificities = c()
thresholds = c()

#par(mfrow=c(3,3))
set.seed(6090)
library(tidyverse)
library(caret)
library(nnet)
library(pROC)
library(caret)
library(randomForest)
library(caret)
library(tree)
library(leaps)
for(r in 1:5){
  n = nrow(pros_new_dat)
  i.cv = sample(1:n,n, replace=FALSE)
  print(i.cv[1])
  pc_copy = pros_new_dat[i.cv,]
  pc.i <- createDataPartition(pc_copy$PathStage,p=0.70, list=F)
  pc_train = pc_copy[pc.i,]
  pc_test = pc_copy[-pc.i,]
  
  
  
  K = 10
  length = length(pc.i)
  folds = cut(1:length, K, labels=FALSE)
  acc.mn = numeric(K)
  feature.matrix = matrix(0,nrow = K, ncol = ncol(pc_train))
  colnames(feature.matrix) = colnames(pc_train)
  feature.matrix = feature.matrix[,-18574]
  
  #consolidated_list = list()
  
  #set.seed(6090)
  for(k in 1:K){
    i.train	= which(folds!=k)
    
    pros_train_old = pc_train[i.train,]
    pros_valid_old = pc_train[-i.train,]
    pc_numeric_pros_train = pros_train_old[,-c(1,2,4,5,6,8,9)]
    pc_categorical_pros_train = pros_train_old[,c(1,2,4,5,6,8,9,18574)]
    
    #chi sqr
    features=c()
    categorical_features=c()
    temp=c()
    adjust=c()
    for(col in 1:(ncol(pc_categorical_pros_train)-1)){
      chisq_test <- chisq.test(pc_categorical_pros_train[,col], pc_categorical_pros_train$PathStage, correct = FALSE)
      temp[col]=chisq_test$p.value
      #if(chisq_test$p.value<0.05){
      #  categorical_features = c(categorical_features,colnames(pc_categorical_pros_train)[col])
      #}
    }
    adjust=p.adjust(temp,method = "fdr")
    for(j in 1:(ncol(pc_categorical_pros_train)-1)){
      if(adjust[j]<0.05){
        categorical_features=c(categorical_features,colnames(pc_categorical_pros_train)[j])
      }
    }
    
    #mutual info
    numerical_features=c()
    mrna.feature.matrix.mmi = matrix(nrow = 1, ncol = ncol(pc_numeric_pros_train)-1)
    colnames(mrna.feature.matrix.mmi) = colnames(pc_numeric_pros_train[1:18566])
    for(i in 1:18566){
      mmi.val=mmi(as.matrix(pc_numeric_pros_train[,i]),data.frame(pc_numeric_pros_train$PathStage))
      mrna.feature.matrix.mmi[1,i]=mmi.val$mi
    }
    mmi.val09.df=data.frame(mrna.feature.matrix.mmi)
    mmi.val09.df=mmi.val09.df[colSums(mmi.val09.df>0.13)>0]
    numerical_features=c(numerical_features,colnames(mmi.val09.df))
    
    
    features=c(numerical_features,categorical_features,"PathStage")
    cat("r = ",r,"k = ",k,"sel features = ",features)
    
    pros_train=pros_train_old[(colnames(pros_train_old) %in% features)]
    pros_valid=pros_valid_old[(colnames(pros_valid_old) %in% features)]
    
    
    reg.fwd = regsubsets(PathStage~., data=pros_train, method="forward")
    R2adj = summary(reg.fwd)$adjr2
    R2adj.index = which.max(R2adj)
    mod = summary(reg.fwd)$which[R2adj.index,]
    feature.matrix[k,names(which(mod)[-1]-1)] = 1
    #feature.matrix[k,-(which(mod)[-1]-1)] = 0
    vars = rownames(data.frame(which(mod)[-1]-1))
    pros_fs_train = pros_train[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_train)))]
    pros_fs_valid = pros_valid[,c(match(vars,colnames(pros_train)),match('PathStage',colnames(pros_valid)))]
    #mn.o = nnet::multinom(PathStage~., data=pros_fs_train)
    mn.o=glm(PathStage~., data=pros_fs_train,family = binomial(logit))
    mn.p = predict(mn.o, newdata=pros_fs_valid,type="response")
    predicted.classes <- ifelse(mn.p < 0.5, "Advanced", "Early")
    tb.mn = table(predicted.classes, pros_fs_valid$PathStage)
    acc.mn[k] = sum(diag(tb.mn)) / sum(tb.mn)
    #consolidated_list <- append(consolidated_list, list(c(vars,acc.mn[k])))
    cat("r = ",r,"k = ",k,"internal acc = ",acc.mn[k])
  }
  #print(consolidated_list)
  accuracy.avg.inner[r] = mean(acc.mn)
  cat("r = ",r,"mean accuracy is :",accuracy.avg.inner[r])
  cat("r = ",r,"SD of accuracy is :",sd(acc.mn))
  fs = c()
  for(i in 1:ncol(feature.matrix)){
    if((sum(feature.matrix[,i])/nrow(feature.matrix)) >= 0.5){
      fs = c(fs, colnames(feature.matrix)[i])
    }
  }
  fs = c(fs, "PathStage")
  print(fs)
  #fs=unique(fs)
  #print(fs)
  pc_train_fs = pc_train[,fs] 
  pc_test_fs = pc_test[,fs] 
  mn.fit = glm(PathStage~., data=pc_train_fs,family = binomial(logit))
  mn.pred = predict(mn.fit, newdata=pc_test_fs,type="response")
  #acc.dt = sum(diag(tb.dt)) / sum(tb.dt)
  #table.mn = confusionMatrix(mn.pred, pc_test_fs$PathStage, mode = 'everything')
  predicted.classes.out <- ifelse(mn.pred < 0.5, "Advanced", "Early")
  tb.mn.out = table(predicted.classes.out, pc_test_fs$PathStage)
  accuracy.outer[r] = round(sum(diag(tb.mn.out)) / sum(tb.mn.out),2)
  #glm.p = predict(mn.fit, pc_test_fs, type="prob")[,2] 
  roc.glm = roc(response=pc_test_fs$PathStage, predictor=mn.pred)
  #roc.pred=as.numeric(predict(mn.fit,pc_test_fs,type='class'))
  #roc.multi=multiclass.roc(pc_test_fs$PathStage,roc.pred)
  #rs=roc.multi[['rocs']]
  #ggpl=ggroc(rs)
  #rs.list=append(rs.list,ggpl)
  precision[r]=round(tb.mn.out[1]/(tb.mn.out[1]+tb.mn.out[3]),2)
  recall[r]=round(tb.mn.out[1]/(tb.mn.out[1]+tb.mn.out[2]),2)
  #precision2=tb.mn.out[1]/(tb.mn.out[1]+tb.mn.out[3])
  #precision[r] = round(table.mn$byClass[5],2)
  #recall[r] = round(table.mn$byClass[6],2)
  f1.score[r] = round(2/((1/precision[r])+(1/recall[r])),2)
  #rf.p = predict(rf.fit, pc_test_fs, type="prob")[,2] 
  ####********need to find precision**************************
  #roc.rf = roc(response=pc_test_fs$PathStage, predictor=rf.p)
  roc.outer[[r]] = roc.glm
  auc.outer[r] = round(roc.glm$auc,2)
  sensitivities[r] = list(roc.glm$sensitivities)
  specificities[r] = list(roc.glm$specificities)
  thresholds[r] = list(roc.glm$thresholds)
}

plot(roc.outer[[1]], col='black')
plot(roc.outer[[2]], add=TRUE, col='red')
plot(roc.outer[[3]], add=TRUE, col='skyblue')
plot(roc.outer[[4]], add=TRUE, col='green')
plot(roc.outer[[5]], add=TRUE, col='orange')
legend(x = 'bottomright', legend=c("Iter 1", "Iter 2", "Iter 3","Iter 4","Iter 5"), 
       col = c("black","red",'skyblue', 'green','orange'), lty=1, cex=1.0, title = 'GLM'
)

max_len_sensitivity = max(lengths(sensitivities))

for(i in 1:length(sensitivities)){
  if(length(sensitivities[[i]])<max_len_sensitivity){
    repeats = max_len_sensitivity-length(sensitivities[[i]])
    sensitivities[[i]] = append(sensitivities[[i]], rep(sensitivities[[i]][length(sensitivities[[i]])],times=repeats))
  }
}

max_len_specificity = max(lengths(specificities))

for(i in 1:length(specificities)){
  if(length(specificities[[i]])<max_len_specificity){
    repeats = max_len_specificity-length(specificities[[i]])
    specificities[[i]] = append(specificities[[i]], rep(specificities[[i]][length(specificities[[i]])],times=repeats))
  }
}


new_sensitivities = c()

for(i in 1:max_len_sensitivity){
  #mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i])/2
  
  mean_sensitivity_i = (sensitivities[[1]][i]+sensitivities[[2]][i]+sensitivities[[3]][i]+sensitivities[[4]][i]+sensitivities[[5]][i])/5
  new_sensitivities = c(new_sensitivities,mean_sensitivity_i)
}
new_sensitivities

new_specificities = c()

for(i in 1:max_len_specificity){
  #mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i])/2
  
  mean_specificity_i = (specificities[[1]][i]+specificities[[2]][i]+specificities[[3]][i]+specificities[[4]][i]+specificities[[5]][i])/5
  new_specificities = c(new_specificities,mean_specificity_i)
}
new_specificities
par(pty='s')
my_labels <- sprintf(rev(seq(0.0, 1.0, 0.2)),
                     fmt = '%#.1f')


par(pty='s')
plot(1-new_specificities, new_sensitivities, type='l',xlab='1-Specificity',ylab='Sensitivity')
abline(coef=c(0,1))





