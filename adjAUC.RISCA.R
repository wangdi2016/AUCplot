library(RISCA)

## 
## The code is for AUC plot (Difei Wang, 2024-05-03)
## 
## the following sentences are from Minta Thomas from Li Hsu group at Fred Hutch.

##
# Please see below the code for calculation age and sex adjusted AUC:
# Input file: a data frame with vectors, age,sex, PRS, outcome and PCs

## Step 1:
## lm1 <- lm(GRS3 ~ ageCRC + sex+pc1+pc2+pc3+pc4 , data=ffdata[ffdata$outcome == 0,]) # you can exclude PCs if you do not use them
#GRS3 is the PRS

## Step 2:
## ffdata$prs_std <- (ffdata$GRS3 - (lm1$coef[1] + lm1$coef[2] * ffdata$ageCRC +
##                                                 lm1$coef[3] * (ffdata$sex)  +
##                                                 lm1$coef[4] * ffdata$pc1    +
##                                                 lm1$coef[5] * ffdata$pc2    +
##                                                 lm1$coef[6] * ffdata$pc3    +
##                                                 lm1$coef[7] * ffdata$pc4 )) / sd(lm1$residuals)

## Step 3:
## adjusted.ROC(status="outcome", variable="prs_std",
## confounders=~ageCRC+sex+pc1+pc2+pc3+pc4 , database=ffdata )$auc

#Please refer ROCt package for the details, alternately you can use RISCA R package with roc.binary() and let us know if you have any questions.


###
### I used RISCA package to generate the AUC plot for CRC project.
### The following code was written with the help of above mention code from Minta and Li.
###

### read table
dat=read.table(file="all.ukbEURu_hm3_shrunk_sparse_mldm.p0.00001.profile.plcoid.entryagebq.sex.age.eu_pca.MAT.clean.txt",sep="\t",header=T)

### PC1-10
dat$pc1 = dat$EU_PCA_component1
dat$pc2 = dat$EU_PCA_component2
dat$pc3 = dat$EU_PCA_component3
dat$pc4 = dat$EU_PCA_component4
dat$pc5 = dat$EU_PCA_component5
dat$pc6 = dat$EU_PCA_component6
dat$pc7 = dat$EU_PCA_component7
dat$pc8 = dat$EU_PCA_component8
dat$pc9 = dat$EU_PCA_component9
dat$pc10 = dat$EU_PCA_component10

###
### crude estimate
###
roc0 <- roc.binary(status="j_colo_cancer", variable="SCORESUM", confounders=~1, data=dat, precision=seq(0.1,0.9, by=0.1))

###
### adjust sex and age
###
lm1 <- lm(SCORESUM ~ age+sex, data=dat) #SCORESUM is the PRS
dat$prs_std <- (dat$SCORESUM - (lm1$coef[1] + lm1$coef[2]  * dat$age +
                                              lm1$coef[3]  * dat$sex))/sd(lm1$residuals) 
roc1 <- roc.binary(status="j_colo_cancer", variable="prs_std", confounders=~age+sex, data=dat, precision=seq(0.1,0.9, by=0.1))

###
### adjust sex, age and 10 EU_PCA
### 
lm2 <- lm(SCORESUM ~ age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data=dat) #SCORESUM is the PRS, 10 eu_pca_component

dat$prs_std_pca <- (dat$SCORESUM - (lm2$coef[1] + lm2$coef[2]  * dat$age + 
                                                  lm2$coef[3]  * dat$sex + 
                                                  lm2$coef[4]  * dat$pc1 +
                                                  lm2$coef[5]  * dat$pc2 +
                                                  lm2$coef[6]  * dat$pc3 +
                                                  lm2$coef[7]  * dat$pc4 + 
                                                  lm2$coef[8]  * dat$pc5 + 
                                                  lm2$coef[9]  * dat$pc6 + 
                                                  lm2$coef[10] * dat$pc7 + 
                                                  lm2$coef[11] * dat$pc8 + 
                                                  lm2$coef[12] * dat$pc9 + 
                                                  lm2$coef[13] * dat$pc10))/sd(lm2$residuals)

roc2 <- roc.binary(status="j_colo_cancer", variable="prs_std_pca", confounders=~age+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data=dat, precision=seq(0.1,0.9, by=0.1))

###
### The corresponding ROC graph
###
pdf(file="col.PRS.adjAUC.RISCA.pdf",height=10, width=20)
par(mfrow=c(1,2))

plot(roc0, type="b", col=1, pch=1, lty=1, ylab="sensitivity", xlab="1-specificity")
lines(roc1, type="b", col=2, pch=2, lty=2)
lines(roc2, type="b", col=3, pch=3, lty=3)

#abline(c(0,0), c(1,1), lty=2, lwd=2)
legend("bottomright", lty=1:3, cex=1.5, lwd=2, col=1:3, pch=1:3, c(paste("no adjust,     (AUC=", round(roc0$auc, 4), ")", sep=""),
                                                                   paste("Adjusted(sex,age),     (AUC=", round(roc1$auc, 4), ")", sep=""),
                                                                   paste("Adjusted(sex,age,Eu-PCA), (AUC=", round(roc2$auc, 4), ")", sep="")))
title("AUC for PLCO Colon Cancer(European) PRS\nall.ukbEURu_hm3_shrunk_sparse_mldm.p0.00001.profile")

### plot density plots
dat1 <- dat[dat$j_colo_cancer==0,]
dat2 <- dat[dat$j_colo_cancer==1,]

dx <- density(dat1$SCORESUM) ## control
dy <- density(dat2$SCORESUM) ## case

plot(dx, col = "blue", lwd = 2, lty=1, ## control
    main = "Density of PLCO European Case and Control PRS", xlab = "",
    xlim = c(min(dx$x, dy$x), c(max(dx$x, dy$x))),  # Min and Max X-axis limits
    ylim = c(min(dx$y, dy$y), c(max(dx$y, dy$y))))  # Min and Max Y-axis limits
lines(dy, col = "red", lwd = 2, lty=2) ## case
legend("topright", c("control", "case"), col = c("blue", "red"), lty = 1:2, cex=1.5, lwd=3)

dev.off()

