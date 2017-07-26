


##################################################################  SIMULATION   #######################################################################################################
rm(list=ls())
setwd("path")
load(".RData")
library(lme4)
library(sva)

totalsamples <- 251 
Simulations <- 100 #the analysis was divided into 10 scripts by changing the set.seed() function below
genesSig <- 500 #true positives
SD.Add <- 0.1 #effect size
set.seed(1)

Matrix.Results_N <- c()
Matrix.Results_mean <- c()
Matrix.Results_pval <- c()


for (k in 1:Simulations){
print(k)
Samples <- c()

coefs<-cbind(rep(NA,nrow(exp.data)), rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)))
colnames(coefs)<-c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
row.names(coefs) <- rownames(exp.data)

coefslm<-cbind(rep(NA,nrow(exp.data)), rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)))
colnames(coefslm)<-c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
row.names(coefslm) <- rownames(exp.data)

coefslmBatch<-cbind(rep(NA,nrow(exp.data)), rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)))
colnames(coefslmBatch)<-c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
row.names(coefslmBatch) <- rownames(exp.data)

coefscom<-cbind(rep(NA,nrow(exp.data)), rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)),rep(NA,nrow(exp.data)))
colnames(coefscom)<-c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
row.names(coefscom) <- rownames(exp.data)


## QCs
 
batchEffect <- rnorm(14,0,gamma)

	QC1.day1 <- QCs.data[,"X445490_2_3_QC1"] #name of the quality control samples in the file
	QC2.day1 <- QCs.data[,"X445490_2_3_QC1"]

	QC1.day2 <- QCs.data[,"X446530_2_3_QC16"]
	QC2.day2 <- QCs.data[,"X446530_2_4_QC17"]

	QC1.day3 <- QCs.data[,"X446521_2_3_QC18"]
	QC2.day3 <- QCs.data[,"X446521_2_4_QC19"]

	QC1.day4 <- QCs.data[,"X445493_2_3_QC3"]
	QC2.day4 <- QCs.data[,"X445493_2_4_QC4"]

	QC1.day5 <- QCs.data[,"X445496_2_3_QC5"]
	QC2.day5 <- QCs.data[,"X445496_2_4_QC6"]

	QC1.day6 <- QCs.data[,"X446524_2_3_QC20"]
	QC2.day6 <- QCs.data[,"X446524_2_4_QC21"]

	QC1.day7 <- QCs.data[,"X445511_2_3_QC7"]
	QC2.day7 <- QCs.data[,"X445511_2_4_QC8"]

 	QC1.day8 <- QCs.data[,"X446511_2_3_QC9"]
	QC2.day8 <- QCs.data[,"X446511_2_4_QC10"] 

	QC1.day9 <- QCs.data[,"X446572_2_3_QC26"]
	QC2.day9 <- QCs.data[,"X446572_2_4_QC27"]

	QC1.day10 <- QCs.data[,"X446514_2_3_QC11"]
	QC2.day10 <- QCs.data[,"X446514_2_4_QC12"]

	QC1.day11 <- QCs.data[,"X446610_2_2_A_QC28"]
	QC2.day11 <- QCs.data[,"X446610_2_3_A_QC29"]

	QC1.day12 <- QCs.data[,"X446517_2_3_QC13"]
	QC2.day12 <- QCs.data[,"X446517_2_4_QC14"]

 	QC1.day13 <- QCs.data[,"X446520_1_3_Qc15"]

 	QC1.day14 <- QCs.data[,"X446577_2_3_A_BaQC31"]
	QC2.day14 <- QCs.data[,"X446577_2_4_A_BaQC32"]


	
#Calculate new QC scaling factor
average.day1.Sim=rowMeans(cbind(QC1.day1,QC2.day1), na.rm=TRUE)
average.day2.Sim=rowMeans(cbind(QC1.day2,QC2.day2), na.rm=TRUE)
average.day3.Sim=rowMeans(cbind(QC1.day3,QC2.day3), na.rm=TRUE)
average.day4.Sim=rowMeans(cbind(QC1.day4,QC2.day4), na.rm=TRUE)
average.day5.Sim=rowMeans(cbind(QC1.day5,QC2.day5), na.rm=TRUE)
average.day6.Sim=rowMeans(cbind(QC1.day6,QC2.day6), na.rm=TRUE)
average.day7.Sim=rowMeans(cbind(QC1.day7,QC2.day7), na.rm=TRUE)
average.day8.Sim=rowMeans(cbind(QC1.day8,QC2.day8), na.rm=TRUE)
average.day9.Sim=rowMeans(cbind(QC1.day9,QC2.day9), na.rm=TRUE)
average.day10.Sim=rowMeans(cbind(QC1.day10,QC2.day10), na.rm=TRUE)
average.day11.Sim=rowMeans(cbind(QC1.day11,QC2.day11), na.rm=TRUE)
average.day12.Sim=rowMeans(cbind(QC1.day12,QC2.day12), na.rm=TRUE)
average.day13.Sim=QC1.day13
average.day14.Sim=rowMeans(cbind(QC1.day14,QC2.day14), na.rm=TRUE)

AllQCs.Sim=rowMeans(cbind(QC1.day1,QC2.day1,QC1.day2,QC2.day2,QC1.day3,QC2.day3,QC1.day4,QC2.day4,QC1.day5,QC2.day5,QC1.day6,QC2.day6,QC1.day7,QC2.day7,QC1.day8,QC2.day8,QC1.day9,QC2.day9,QC1.day10,QC2.day10,QC1.day11,QC2.day11,QC1.day12,QC2.day12,QC1.day13,QC1.day14,QC2.day14), na.rm=TRUE)

factor1.Sim <-average.day1.Sim-AllQCs.Sim
factor2.Sim <-average.day2.Sim-AllQCs.Sim
factor3.Sim <-average.day3.Sim-AllQCs.Sim
factor4.Sim <-average.day4.Sim-AllQCs.Sim
factor5.Sim <-average.day5.Sim-AllQCs.Sim
factor6.Sim <-average.day6.Sim-AllQCs.Sim
factor7.Sim <-average.day7.Sim-AllQCs.Sim
factor8.Sim <-average.day8.Sim-AllQCs.Sim
factor9.Sim <-average.day9.Sim-AllQCs.Sim
factor10.Sim <-average.day10.Sim-AllQCs.Sim
factor11.Sim <-average.day11.Sim-AllQCs.Sim
factor12.Sim <-average.day12.Sim-AllQCs.Sim
factor13.Sim <-average.day13.Sim-AllQCs.Sim
factor14.Sim <-average.day14.Sim-AllQCs.Sim



Add.subject <- c()
Add.metadata <- c()
for (m in 1:totalsamples){ 
Add <- rnorm(1,0,SD.Add)
u <- c(rep(Add,genesSig),rep(0,nrow(exp.data)-genesSig))
Add.subject <- cbind(Add.subject,u)
Add.metadata <- c(Add.metadata,Add)

}

condition.day1 <- c()
condition.day2 <- c()
condition.day3 <- c()
condition.day4 <- c()
condition.day5 <- c()
condition.day6 <- c()
condition.day7 <- c()
condition.day8 <- c()
condition.day9 <- c()
condition.day10 <- c()
condition.day11 <- c()
condition.day12 <- c()
condition.day13 <- c()
condition.day14 <- c()


	for (j in 2:21){ 

			SX.dayX <- exp.data[,j] + factor1.Sim  + rep(batchEffect[1],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day1 <- names(exp.data[j])
			condition.day1 <- rbind(condition.day1,names.day1)
	}

	for (j in 22:43){ 

			SX.dayX <- exp.data[,j] + factor2.Sim  + rep(batchEffect[2],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day2 <- names(exp.data[j])
			condition.day2 <- rbind(condition.day2,names.day2)
	}

	for (j in 44:64){ 

			SX.dayX <- exp.data[,j] + factor3.Sim  + rep(batchEffect[3],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day3 <- names(exp.data[j])
			condition.day3 <- rbind(condition.day3,names.day3)
	}

	for (j in 65:85){ 

			SX.dayX <- exp.data[,j] + factor4.Sim  + rep(batchEffect[4],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day4 <- names(exp.data[j])
			condition.day4 <- rbind(condition.day4,names.day4)
	}

	for (j in 86:106){ 

			SX.dayX <- exp.data[,j] + factor5.Sim  + rep(batchEffect[5],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day5 <- names(exp.data[j])
			condition.day5 <- rbind(condition.day5,names.day5)
	}

	for (j in 107:127){ 

			SX.dayX <- exp.data[,j] + factor6.Sim  + rep(batchEffect[6],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day6 <- names(exp.data[j])
			condition.day6 <- rbind(condition.day6,names.day6)
	}

	for (j in 128:147){ 

			SX.dayX <- exp.data[,j] + factor7.Sim  + rep(batchEffect[7],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day7 <- names(exp.data[j])
			condition.day7 <- rbind(condition.day7,names.day7)
	}

	for (j in 167:186){ 

			SX.dayX <- exp.data[,j] + factor8.Sim  + rep(batchEffect[8],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day8 <- names(exp.data[j])
			condition.day8 <- rbind(condition.day8,names.day8)
	}

	for (j in 207:228){ 

			SX.dayX <- exp.data[,j] + factor9.Sim  + rep(batchEffect[9],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day9 <- names(exp.data[j])
			condition.day9 <- rbind(condition.day9,names.day9)
	}

	for (j in 229:248){ 

			SX.dayX <- exp.data[,j] + factor10.Sim  + rep(batchEffect[10],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day10 <- names(exp.data[j])
			condition.day10 <- rbind(condition.day10,names.day10)
	}

	for (j in 249:260){ 

			SX.dayX <- exp.data[,j] + factor11.Sim  + rep(batchEffect[11],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day11 <- names(exp.data[j])
			condition.day11 <- rbind(condition.day11,names.day11)
	}

	for (j in 261:281){ 

			SX.dayX <- exp.data[,j] + factor12.Sim  + rep(batchEffect[12],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day12 <- names(exp.data[j])
			condition.day12 <- rbind(condition.day12,names.day12)
	}

	for (j in 282:284){ 

			SX.dayX <- exp.data[,j] + factor13.Sim  + rep(batchEffect[13],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day13 <- names(exp.data[j])
			condition.day13 <- rbind(condition.day13,names.day13)
	}

	for (j in 285:291){ 

			SX.dayX <- exp.data[,j] + factor14.Sim  + rep(batchEffect[14],nrow(exp.data)) + Add.subject[,j]
			Samples <- cbind(Samples,SX.dayX)
			names.day14 <- names(exp.data[j])
			condition.day14 <- rbind(condition.day14,names.day14)
	}


condition <- rbind(condition.day1,condition.day2,condition.day3,condition.day4,condition.day5,condition.day6,condition.day7,condition.day8,condition.day9,condition.day10,condition.day11,condition.day12,condition.day13, condition.day14)
Add.metadata <- Add.metadata[c(-1,-c(148:166),-c(187:206))]
condition <- cbind.data.frame(condition,Add.metadata,rownames(condition))
colnames(condition) <- c("ID","treatment","batch")
Samples_t <- t(Samples)


##### STATISTICAL MODELS


modcombat = model.matrix(~1, data=condition)
combat_edata = ComBat(Samples, batch=condition$batch, mod=modcombat, par.prior=TRUE)

	for(i in 1:ncol(Samples_t)){

		  LMM <- lmer(Samples_t[,i] ~ treatment + (1|batch), data=condition)
		  out <- coef(summary(LMM))
		  out2 <- cbind(out,`Pr(>|t|)`=2 * pt(abs(out[,"t value"]), df=df.residual(LMM),lower.tail=FALSE))
		  coefs[i,1]<- out2[2,1]
		  coefs[i,2]<- out2[2,2]
		  coefs[i,3]<- out2[2,3]
		  coefs[i,4]<- out2[2,4]
		  
	          LM <- lm(Samples_t[,i] ~ treatment, data=condition)
		  outlm <- coef(summary(LM))
		  out2lm <- cbind(outlm,`Pr(>|t|)`=2 * pt(abs(outlm[,"t value"]), df=df.residual(LM),lower.tail=FALSE))
	          coefslm[i,1]<- out2lm[2,1]
		  coefslm[i,2]<- out2lm[2,2]
		  coefslm[i,3]<- out2lm[2,3]
		  coefslm[i,4]<- out2lm[2,4]

	          LMBatch <- lm(Samples_t[,i] ~ treatment + batch, data=condition)
		  outlmBatch <- coef(summary(LMBatch))
		  out2lmBatch <- cbind(outlmBatch,`Pr(>|t|)`=2 * pt(abs(outlmBatch[,"t value"]), df=df.residual(LMBatch),lower.tail=FALSE))
	          coefslmBatch[i,1]<- out2lmBatch[2,1]
		  coefslmBatch[i,2]<- out2lmBatch[2,2]
		  coefslmBatch[i,3]<- out2lmBatch[2,3]
		  coefslmBatch[i,4]<- out2lmBatch[2,4]
		
		  LMcom <- lm(t(combat_edata)[,i] ~ treatment, data=condition)
		  outcom <- coef(summary(LMcom))
		  out2com <- cbind(outcom,`Pr(>|t|)`=2 * pt(abs(outcom[,"t value"]), df=df.residual(LMcom),lower.tail=FALSE))
	          coefscom[i,1]<- out2com[2,1]
		  coefscom[i,2]<- out2com[2,2]
		  coefscom[i,3]<- out2com[2,3]
		  coefscom[i,4]<- out2com[2,4]
	}

pval <- coefs[,"Pr(>|t|)"]
pval_fdr<-p.adjust(pval, method="fdr")
pval_bonf<-p.adjust(pval, method="bonferroni")

pvallm <- coefslm[,"Pr(>|t|)"]
pval_fdrlm<-p.adjust(pvallm, method="fdr")
pval_bonflm<-p.adjust(pvallm, method="bonferroni")

pvallmBatch <- coefslmBatch[,"Pr(>|t|)"]
pval_fdrlmBatch<-p.adjust(pvallmBatch, method="fdr")
pval_bonflmBatch<-p.adjust(pvallmBatch, method="bonferroni")

pvalcom <- coefscom[,"Pr(>|t|)"]
pval_fdrcom<-p.adjust(pvalcom, method="fdr")
pval_bonfcom<-p.adjust(pvalcom, method="bonferroni")


TP_LMM <- length(which (pval_fdr[1:genesSig] <0.05))  
FP_LMM <- length(which (pval_fdr[(genesSig+1):nrow(exp.data)] <0.05))

TP_LM <- length(which (pval_fdrlm[1:genesSig] <0.05))
FP_LM <- length(which (pval_fdrlm[(genesSig+1):nrow(exp.data)] <0.05))

TP_LMBatch <- length(which (pval_fdrlmBatch[1:genesSig] <0.05))
FP_LMBatch <- length(which (pval_fdrlmBatch[(genesSig+1):nrow(exp.data)] <0.05))

TP_LMcom <- length(which (pval_fdrcom[1:genesSig] <0.05))
FP_LMcom <- length(which (pval_fdrcom[(genesSig+1):nrow(exp.data)] <0.05))

Results_N <- cbind(TP_LMM,FP_LMM,TP_LM,FP_LM,TP_LMBatch,FP_LMBatch,TP_LMcom,FP_LMcom) 
Matrix.Results_N <- rbind(Matrix.Results_N,Results_N)

TP_LMM_mean <- mean(pval_fdr[1:genesSig])  
FP_LMM_mean <- mean(pval_fdr[(genesSig+1):nrow(exp.data)] )

TP_LM_mean <- mean(pval_fdrlm[1:genesSig] )
FP_LM_mean<- mean(pval_fdrlm[(genesSig+1):nrow(exp.data)] )

TP_LMBatch_mean <- mean(pvallmBatch[1:genesSig] )
FP_LMBatch_mean <- mean(pvallmBatch[(genesSig+1):nrow(exp.data)] )

TP_LMcom_mean <- mean(pval_fdrcom[1:genesSig] )
FP_LMcom_mean <- mean(pval_fdrcom[(genesSig+1):nrow(exp.data)] )

Results_mean <- cbind(TP_LMM_mean,FP_LMM_mean,TP_LM_mean,FP_LM_mean,TP_LMBatch_mean,FP_LMBatch_mean,TP_LMcom_mean,FP_LMcom_mean) 
Matrix.Results_mean <- rbind(Matrix.Results_mean,Results_mean)


TP_LMM_mean_pval <- mean(pval[1:genesSig])  
FP_LMM_mean_pval <- mean(pval[(genesSig+1):nrow(exp.data)] )

TP_LM_mean_pval <- mean(pvallm[1:genesSig] )
FP_LM_mean_pval<- mean(pvallm[(genesSig+1):nrow(exp.data)] )

TP_LMBatch_mean_pval <- mean(pvallmBatch[1:genesSig] )
FP_LMBatch_mean_pval <- mean(pvallmBatch[(genesSig+1):nrow(exp.data)] )

TP_LMcom_mean_pval <- mean(pval_fdrcom[1:genesSig] )
FP_LMcom_mean_pval <- mean(pval_fdrcom[(genesSig+1):nrow(exp.data)] )

Results_mean_pval <- cbind(TP_LMM_mean_pval,FP_LMM_mean_pval,TP_LM_mean_pval,FP_LM_mean_pval,TP_LMBatch_mean_pval,FP_LMBatch_mean_pval,TP_LMcom_mean_pval,FP_LMcom_mean_pval) 
Matrix.Results_pval <- rbind(Matrix.Results_pval,Results_mean_pval)

}

write.table(Matrix.Results_N,paste("Results Number FP-TP 1. Runnings:",k,".txt",sep=''), sep="\t", quote=F,row.names=TRUE, col.names=NA)
write.table(Matrix.Results_mean,paste("Results mean FP-TP 1. Runnings:",k,".txt",sep=''), sep="\t", quote=F,row.names=TRUE, col.names=NA)
write.table(Matrix.Results_pval,paste("Results mean pval FP-TP 1. Runnings:",k,".txt",sep=''), sep="\t", quote=F,row.names=TRUE, col.names=NA)



