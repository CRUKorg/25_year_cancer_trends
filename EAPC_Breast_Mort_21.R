##*************************************************************************
##*************************************************************************
##**                                                              	***
##**		Cancer Intelligence - CRUK						***
##**												***
##**		Filename: EAPC_Breast_Mort.R						***
##**												***
##**		Purpose of prog.:								***
##**		calculate the Estimated Annual Percentage Change (EAPC)	***
##**		of rates during a time period with the Confidence		***
##**		Interval (CI)				 				***
##** 												***
##**		Operating System: MS Windows/ R					***
##**		Version: 1/ A. S. Ahmad / 03.06.2019	             	***
##**												***
##*************************************************************************
##*************************************************************************

##** Clear all
  rm(list=ls(all=TRUE))

##** Get data
  Dat35 <- read.csv("G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\Raw\\ASR_BreastMort_35to49_to_2018.csv")
  Dat50 <- read.csv("G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\Raw\\ASR_BreastMort_50to69_to_2018.csv")
  Dat <- rbind(Dat35,Dat50)
  rm(Dat35);rm(Dat50)


##** Keep Female data
  vL <- c("Year","Female","AgeGroup")
  dat <- Dat[,vL]
  rm(Dat)

 colnames(dat) <- c("year","asr","age")

 rm(vL)


##** Call library for 
  library(Rcan)
  library(rms)


APC.LR <- APC.QR <- APC.PR <- APC.QPR  <- matrix(NA,nrow=length(unique(dat$age)),ncol=3)

colnames(APC.LR) <- c("eapc","eapc_low","eapc_up")
colnames(APC.QR) <- colnames(APC.PR) <- colnames(APC.QPR) <- colnames(APC.LR)

rownames(APC.QR) <- c("QR_Age:35to49","QR_Age:50to69")
rownames(APC.PR) <- c("PR_Age:35to49","PR_Age:50to69")
rownames(APC.QPR) <- c("QPR_Age:35to49","QPR_Age:50to69")
rownames(APC.LR) <- c("LR_Age:35to49","LR_Age:50to69")

for(i in 1:length(unique(dat$age))){
 fitLR <- glm(log(asr)~year,subset=age==unique(age)[i],data=dat)
 APC.LR[i,] <- c((exp(fitLR$coef[2])-1)*100,
                 (exp(fitLR$coef[2]-qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100,
                 (exp(fitLR$coef[2]+qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100)
fitQR <- Rq(log(asr)~year,subset=age==unique(age)[i],data=dat)
APC.QR[i,] <- c((exp(fitQR$coef[2])-1)*100,
                (exp(fitQR$coef[2]-qnorm(0.975)*sqrt(diag(fitQR$var))[2])-1)*100,
                (exp(fitQR$coef[2]+qnorm(0.975)*sqrt(diag(fitQR$var))[2])-1)*100)

fitPR <- glm(asr~year,family=poisson,subset=age==unique(age)[i],data=dat)
APC.PR[i,] <- c((exp(fitPR$coef[2])-1)*100,
                (exp(fitPR$coef[2]-qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100,
                (exp(fitPR$coef[2]+qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100)

fitQPR <- glm(asr~year,family=quasipoisson,subset=age==unique(age)[i],data=dat)
APC.QPR[i,] <- c((exp(fitQPR$coef[2])-1)*100,
                (exp(fitQPR$coef[2]-qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100,
                (exp(fitQPR$coef[2]+qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100)

}

 APC.Res <- rbind(APC.LR,APC.QR,APC.PR,APC.QPR)

 Mod <- csu_eapc(dat,"asr", "year",group_by="age")




##** Define path where results should be saved
  Path <- "G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\R outputs\\"


write.csv(Mod,paste(Path,paste("Rcan_ASR_BreastMort","csv",sep="."),sep=""),row.names=F)
write.csv(APC.Res,paste(Path,paste("Res_ASR_BreastMort","csv",sep="."),sep=""))
