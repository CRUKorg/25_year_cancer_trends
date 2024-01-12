##*************************************************************************
##*************************************************************************
##**                                                              	***
##**		Cancer Intelligence - CRUK						***
##**												***
##**		Filename: EAPC_Female.R							***
##**												***
##**		Purpose of prog.:								***
##**		calculate the Estimated Annual Percentage Change (EAPC)	***
##**		of rates during a time period with the Confidence		***
##**		Interval (CI)				 				***
##** 												***
##**		Operating System: MS Windows/ R					***
##**		Version: 1/ A. S. Ahmad / 30.05.2019	             	***
##**												***
##*************************************************************************
##*************************************************************************

library("dplyr")
library("tidyverse")

##** Clear all
  rm(list=ls(all=TRUE))

##** Get data
  Dat <- read.csv("G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\Raw\\Table3&4 inc_mort_93to18.csv",  fileEncoding = 'UTF-8-BOM')

##** Keep Female data
  vL <- c("CancerSiteLevelThreeCodeAndName","ThreeYearPeriod",
         "IncidenceASRFemaleThreeYrRolling","MortalityASRFemaleThreeYrRolling")
  
  datx <- select(Dat, c("CancerSiteLevelThreeCodeAndName","ThreeYearPeriod",
                        "IncidenceASRFemaleThreeYrRolling","MortalityASRFemaleThreeYrRolling"))
  rm(Dat)

##** AS characters
  datx$ThreeYearPeriod <- as.character(datx$ThreeYearPeriod)
  datx$CancerSiteLevelThreeCodeAndName <- as.character(datx$CancerSiteLevelThreeCodeAndName)

##** remove C62 and C62 from the statistical analysis
  dat <- datx[!datx$CancerSiteLevelThreeCodeAndName%in%c("C62: Testis","C61: Prostate"),]
  dim(dat)
  dim(datx)
 (dim(datx)[1]-dim(dat)[1])/length(1993:2016)
  rm(datx)

##** Get cancer sites
  Site <- unique(dat$CancerSiteLevelThreeCodeAndName)
##** Get type of data ASR Incidence and Mortality
  Type <- c("IncidenceASRFemaleThreeYrRolling","MortalityASRFemaleThreeYrRolling")

##** make year numeric
  yx <- c()
  for(j in 1:dim(dat)[1]){
    yx[j] <- as.numeric(strsplit(dat$ThreeYearPeriod[j],split="-")[[1]][1])
  }
  dat$year <- yx; rm(yx)
  table(dat$year)
  length(Site)

##** Define path where results should be saved
  Path <- "G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\R outputs\\"

##** Call library for
library(Rcan)
library(rms)


##** Loop for the analysis
for(k in 1:2){
dd <- dat[,c("CancerSiteLevelThreeCodeAndName",Type[k],"year")]
colnames(dd) <- c("CancerSiteLevelThreeCodeAndName","asr","year")


Res <- matrix(NA,nrow=length(Site),ncol=3)
rownames(Res) <- Site
colnames(Res) <- c("eapc","eapc_low","eapc_up")


APC.LR <- APC.QR <- APC.PR <- APC.QPR  <- matrix(NA,nrow=length(Site),ncol=3)

rownames(APC.LR) <- Site
colnames(APC.LR) <- c("eapc","eapc_low","eapc_up")

colnames(APC.QR) <- colnames(APC.PR) <- colnames(APC.QPR) <- colnames(APC.LR)
rownames(APC.QR) <- rownames(APC.PR) <- rownames(APC.QPR) <- rownames(APC.LR)


for(i in 1:length(Site)){
dx <- dd[dd$CancerSiteLevelThreeCodeAndName%in%Site[i],]

fitLR <- glm(log(asr)~year,data=dx)
APC.LR[i,] <- c((exp(fitLR$coef[2])-1)*100,
                (exp(fitLR$coef[2]-qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100,
                (exp(fitLR$coef[2]+qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100)

fitQR <- Rq(log(asr)~year,data=dx)
APC.QR[i,] <- c((exp(fitQR$coef[2])-1)*100,
                (exp(fitQR$coef[2]-qnorm(0.975)*sqrt(diag(fitQR$var))[2])-1)*100,
                (exp(fitQR$coef[2]+qnorm(0.975)*sqrt(diag(fitQR$var))[2])-1)*100)

fitPR <- glm(asr~year,family=poisson,data=dx)
APC.PR[i,] <- c((exp(fitPR$coef[2])-1)*100,
                (exp(fitPR$coef[2]-qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100,
                (exp(fitPR$coef[2]+qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100)

fitQPR <- glm(asr~year,family=quasipoisson,data=dx)
APC.QPR[i,] <- c((exp(fitQPR$coef[2])-1)*100,
                (exp(fitQPR$coef[2]-qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100,
                (exp(fitQPR$coef[2]+qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100)

 Fit <- csu_eapc(dx,"asr", "year")
  Res[i,] <- c(Fit$ eapc, Fit$ eapc_low, Fit$ eapc_up)
  rm(Fit)
  rm(dx)
}

write.csv(Res,paste(Path,paste(paste("Rcan",colnames(dat)[(2+k)],sep="_"),"csv",sep="."),sep=""))

write.csv(APC.LR,paste(Path,paste(paste("LR",colnames(dat)[(2+k)],sep="_"),"csv",sep="."),sep=""))
write.csv(APC.QR,paste(Path,paste(paste("QR",colnames(dat)[(2+k)],sep="_"),"csv",sep="."),sep=""))
write.csv(APC.PR,paste(Path,paste(paste("PR",colnames(dat)[(2+k)],sep="_"),"csv",sep="."),sep=""))
write.csv(APC.QPR,paste(Path,paste(paste("QPR",colnames(dat)[(2+k)],sep="_"),"csv",sep="."),sep=""))



##** Clear all
  rm(list=c("APC.LR","APC.QR","APC.PR","APC.QPR","Res"))

cat("Type ->", Type[k],"\n")
flush.console()

}


