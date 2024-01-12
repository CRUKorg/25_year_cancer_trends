##*************************************************************************
##*************************************************************************
##**                                                              	***
##**		Cancer Intelligence - CRUK						***
##**												***
##**		Filename: EAPC_AllCancerExcC44.R					***
##**												***
##**		Purpose of prog.:								***
##**		calculate the Estimated Annual Percentage Change (EAPC)	***
##**		of rates during a time period with the Confidence		***
##**		Interval (CI)				 				***
##** 												***
##**		Operating System: MS Windows/ R					***
##**		Version: 1/ A. S. Ahmad / 19.06.2019	             	***
##**												***
##*************************************************************************
##*************************************************************************

##** Clear all
  rm(list=ls(all=TRUE))

##** Get data
  DatX <- read.csv("G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\Raw\\Table3&4 inc_mort_93to18.csv",  fileEncoding = 'UTF-8-BOM')
  colnames(DatX) <- c("CancerSiteName","IncidenceASRFemale","IncidenceASRMale","Year",                     
                     "CancerSiteName","MortalityASRFemale","MortalityASRMale","X3Year")         
  DatX$CancerSiteName <- as.character(DatX$CancerSiteName)
  DatX$Year <- as.character(DatX$Year)

##** Get cancer sites
  Site <- unique(DatX$CancerSiteName)

##** make year numeric
  yx <- c()
  for(j in 1:dim(DatX)[1]){
    yx[j] <- as.numeric(strsplit(DatX$Year[j],split="-")[[1]][1])
  }
  DatX$year <- yx; rm(yx)
  table(DatX$year)
  

  wL <- c("CancerSiteName","year","IncidenceASRFemale","IncidenceASRMale","MortalityASRFemale","MortalityASRMale")
  Dat <- DatX[DatX$CancerSiteName%in%c("C54-C55: Uterus"),wL]
  rm(wL);rm(DatX)

##** Keep Female data
  vL <- c("MortalityASRFemale","IncidenceASRFemale") 
  
   APC.LR <- APC.QR <- APC.PR <- APC.QPR  <- matrix(NA,nrow=length(vL),ncol=3)

 rownames(APC.LR) <- vL
 colnames(APC.LR) <- c("eapc","eapc_low","eapc_up")

 colnames(APC.QR) <- colnames(APC.PR) <- colnames(APC.QPR) <- colnames(APC.LR)
 rownames(APC.QR) <- rownames(APC.PR) <- rownames(APC.QPR) <- rownames(APC.LR)

 library(rms)

  for(i in 1:length(vL)){
    dx <- Dat[,c("year",vL[i])]
    colnames(dx) <- c("year","asr")
  
  fitLR <- glm(log(asr)~year,data=dx)
  APC.LR[i,] <- c((exp(fitLR$coef[2])-1)*100,
                  (exp(fitLR$coef[2]-qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100,
                  (exp(fitLR$coef[2]+qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100)

  fitPR <- glm(asr~year,family=poisson,data=dx)
  APC.PR[i,] <- c((exp(fitPR$coef[2])-1)*100,
                  (exp(fitPR$coef[2]-qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100,
                  (exp(fitPR$coef[2]+qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100)

  fitQPR <- glm(asr~year,family=quasipoisson,data=dx)
  APC.QPR[i,] <- c((exp(fitQPR$coef[2])-1)*100,
                   (exp(fitQPR$coef[2]-qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100,
                  (exp(fitQPR$coef[2]+qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100)
}

 Res <- rbind(APC.LR,APC.PR,APC.QPR)
 rownames(Res) <- c("MortalityASRFemale_LR","IncidenceASRFemale_LR",
                    "MortalityASRFemale_PR","IncidenceASRFemale_PR",
                    "MortalityASRFemale_QPR","IncidenceASRFemale_QPR")

 write.csv(Res,"G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\R outputs\\ResUterus.csv")


#################################

##** Clear all
  rm(list=ls(all=TRUE))

##** Get data
  DatX <- read.csv("G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\Raw\\Table3&4 inc_mort_93to18.csv",  fileEncoding = 'UTF-8-BOM')
  colnames(DatX) <- c("CancerSiteName","IncidenceASRFemale","IncidenceASRMale","Year",                     
                     "CancerSiteName","MortalityASRFemale","MortalityASRMale","X3Year")         
  DatX$CancerSiteName <- as.character(DatX$CancerSiteName)
  DatX$Year <- as.character(DatX$Year)

##** Get cancer sites
  Site <- unique(DatX$CancerSiteName)

##** make year numeric
  yx <- c()
  for(j in 1:dim(DatX)[1]){
    yx[j] <- as.numeric(strsplit(DatX$Year[j],split="-")[[1]][1])
  }
  DatX$year <- yx; rm(yx)
  table(DatX$year)
  

  wL <- c("CancerSiteName","year","IncidenceASRFemale","IncidenceASRMale","MortalityASRFemale","MortalityASRMale")
  Dat <- DatX[DatX$CancerSiteName%in%c("C00-C97: All Cancers excl C44 and C61 in men and C50 in women"),wL]
  rm(wL);rm(DatX)

##** Keep male data
  vL <- c("IncidenceASRMale","MortalityASRMale","IncidenceASRFemale","MortalityASRFemale")
 
  
   APC.LR <- APC.QR <- APC.PR <- APC.QPR  <- matrix(NA,nrow=length(vL),ncol=3)

 rownames(APC.LR) <- vL
 colnames(APC.LR) <- c("eapc","eapc_low","eapc_up")

 colnames(APC.QR) <- colnames(APC.PR) <- colnames(APC.QPR) <- colnames(APC.LR)
 rownames(APC.QR) <- rownames(APC.PR) <- rownames(APC.QPR) <- rownames(APC.LR)

 library(rms)

  for(i in 1:length(vL)){
    dx <- Dat[,c("year",vL[i])]
    colnames(dx) <- c("year","asr")
  
  fitLR <- glm(log(asr)~year,data=dx)
  APC.LR[i,] <- c((exp(fitLR$coef[2])-1)*100,
                  (exp(fitLR$coef[2]-qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100,
                  (exp(fitLR$coef[2]+qnorm(0.975)*summary(fitLR)$coef[2,2])-1)*100)

  fitPR <- glm(asr~year,family=poisson,data=dx)
  APC.PR[i,] <- c((exp(fitPR$coef[2])-1)*100,
                  (exp(fitPR$coef[2]-qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100,
                  (exp(fitPR$coef[2]+qnorm(0.975)*summary(fitPR)$coef[2,2])-1)*100)

  fitQPR <- glm(asr~year,family=quasipoisson,data=dx)
  APC.QPR[i,] <- c((exp(fitQPR$coef[2])-1)*100,
                   (exp(fitQPR$coef[2]-qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100,
                  (exp(fitQPR$coef[2]+qnorm(0.975)*summary(fitQPR)$coef[2,2])-1)*100)
}

 Res <- rbind(APC.LR,APC.PR,APC.QPR)
 rownames(Res) <- c("IncidenceASRMale_LR","MortalityASRMale_LR","IncidenceASRFemale_LR","MortalityASRFemale_LR",
                    "IncidenceASRMale_PR","MortalityASRMale_PR","IncidenceASRFemale_PR","MortalityASRFemale_PR",
                    "IncidenceASRMale_QPR","MortalityASRMale_QPR","IncidenceASRFemale_QPR","MortalityASRFemale_QPR")

 write.csv(Res,"G:\\Cancer Information\\Statistical Information Team\\Data Analysis and Statistical Info\\Projects\\Working age scoping for JS\\October 2021 Update\\R outputs\\ResAllCancersExclC44andC61inMenAndC50inWomen.csv")


library(Rcan)
dd <- Dat[,c("year",vL[1])]
colnames(dd)[2] <-"asr"
Fit <- csu_eapc(dd,"asr", "year")
