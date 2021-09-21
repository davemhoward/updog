###########
#Genomic Profile Risk Score analysis considering covariates
#Based on Hong Lee & Naomi Wray December 2013 updated January 2014
###########
h2l_R2N <- function(k, r2n, p) {
  # k baseline disease risk
  # r2n Nagelkerkes attributable to genomic profile risk score
  # proportion of sample that are cases
  # calculates proportion of variance explained on the liability scale
  #from ABC at http://www.complextraitgenomics.com/software/
  #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  h2l_R2N <- cc * e * r2n / (1 + cc * e * theta * r2n)
}
actual <- read.table("scoreoutput_genomewide.txt",header=T,sep="") ## CURRENTLY NO HEADER BUT WILL LIKELY CHANGE
phenotype<-read.table("BroadDepressionPhenotypeforMeta.txt",header=F,sep="")
combi<-merge(actual,phenotype,by.x="ID",by.y="V1")

pcs<-readRDS("../updog/ukb_sqc_qc_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.covars.rds")

ri<-merge(combi[which(combi$V3!=-9),],pcs,by.x="ID",by.y="f.eid")
colnames(ri)<-toupper(colnames(ri))


ri$PHENO<-ri$V3
ri$PHENO1<-ri$PHENO-1  ## -1 minus for phenotype    
P<-sum(ri$PHENO1)/length(ri$PHENO1) # proportion of target sample that are cases
N<-length(ri$PHENO1)
NCA<-sum(ri$PHENO1==1)
NCO<-sum(ri$PHENO1==0)

K <- 0.353 ## broad depression across whole UKB

for (i in 4:14) {

ri$SCORE<-ri[,i]
        
### normalize the score
ri$NSCORE<-(ri$SCORE-mean(ri$SCORE))/sd(ri$SCORE)

tstF <- glm(PHENO1 ~ NSCORE + GENOTYPING.ARRAY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = ri,family = binomial(logit)) # logit model
tstR <- glm(PHENO1 ~ GENOTYPING.ARRAY + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = ri,family = binomial(logit)) # logit model
tst0 <- glm(PHENO1 ~ 1 , data = ri,family = binomial(logit)) # logit model
LLF<-logLik(tstF)
LLR<-logLik(tstR)
LL0<-logLik(tst0)
    
CSv<-1-exp((2/N)*(LLR[1]-LLF[1]))
CS<-1-exp((2/N)*(LL0[1]-LLF[1]))
NK0<-CS/(1-exp((2/N)*LL0[1]))
NKv<-CSv/(1-exp((2/N)*LLR[1]))
    
#pvalue
devdiff<-tstR$deviance-tstF$deviance
df<-tstR$df.residual-tstF$df.residual
pval<-pchisq(devdiff,df,lower.tail=F)
    
h2l_r2n <- h2l_R2N(K,NKv,P) #nagelkerkes

print(data.frame(colnames(ri[i]),N,P,NKv,pval,K,h2l_r2n,NCA,NCO))

}
