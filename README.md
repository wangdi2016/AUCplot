# AUCplot

The code is for AUC plot (DW, 2024-05-03) for the CRC PRS project.

the following sentences are from Minta Thomas from Li Hsu group at Fred Hutch.

Please see below the code for calculation age and sex adjusted AUC:
Input file: a data frame with vectors, age,sex, PRS, outcome and PCs

## Step 1:
lm1 <- lm(GRS3 ~ ageCRC + sex+pc1+pc2+pc3+pc4 , data=ffdata[ffdata$outcome == 0,]) # you can exclude PCs if you do not use them
#GRS3 is the PRS

## Step 2:
ffdata$prs_std <- (ffdata$GRS3 - (lm1$coef[1] + lm1$coef[2] * ffdata$ageCRC +
                                                 lm1$coef[3] * (ffdata$sex)  +
                                                 lm1$coef[4] * ffdata$pc1    +
                                                 lm1$coef[5] * ffdata$pc2    +
                                                 lm1$coef[6] * ffdata$pc3    +
                                                 lm1$coef[7] * ffdata$pc4 )) / sd(lm1$residuals)

## Step 3:
adjusted.ROC(status="outcome", variable="prs_std",
confounders=~ageCRC+sex+pc1+pc2+pc3+pc4 , database=ffdata )$auc

#Please refer ROCt package for the details, alternately you can use RISCA R package with roc.binary() and let us know if you have any questions.

