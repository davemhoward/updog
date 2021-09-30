if (!require("BEDMatrix")) {
  install.packages("BEDMatrix", repos = "https://cloud.r-project.org")
} 
library("BEDMatrix")

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
chunk <- as.numeric(args[2])
name <- args[3]

## Pull in chunked test data files
geno<-BEDMatrix(paste0("temp/testdata_",name,"_chr",chr,"_",sprintf("%03d",chunk)))
bim<-read.table(paste0("temp/testdata_",name,"_chr",chr,"_",sprintf("%03d",chunk),".bim"),header=F,sep="")
fam<-read.table(paste0("temp/testdata_",name,"_chr",chr,"_",sprintf("%03d",chunk),".fam"),header=F,sep="")

## Pull in chunked ld files
ld<-BEDMatrix(paste0("temp/ldref_",name,"_chr",chr,"_",sprintf("%03d",chunk)))
ldbim<-read.table(paste0("temp/ldref_",name,"_chr",chr,"_",sprintf("%03d",chunk),".bim"),header=F,sep="",stringsAsFactors=F)
ldbim$ld<-0
window<-250000 ## set window to 250kb

## Pull in chunked summary statistics
sumstats<-read.table(paste0("temp/sumstats_",name,"_chr",chr,"_",sprintf("%03d",chunk),".txt"),header=F,sep="",stringsAsFactors=F)

## Pull in chunked genetic scores
scores<-read.table(paste0("temp/chunkscores_",name,"_chr",chr,"_",sprintf("%03d",chunk)),header=F,sep="",stringsAsFactors=F)

ORIGINALSUMSCORE<-0
UPDOGSUMSCORE1a<-0
UPDOGSUMSCORE1b<-0
UPDOGSUMSCORE1c<-0
UPDOGSUMSCORE1d<-0
UPDOGSUMSCORE1e<-0
riskscoredown<-structure(1:(nrow(geno)), names=rownames(geno))
riskscoreup<-structure(1:(nrow(geno)), names=rownames(geno))
UPDOGSUMSCORE2a<-0
UPDOGSUMSCORE2b<-0
UPDOGSUMSCORE2c<-0
UPDOGSUMSCORE2d<-0
UPDOGSUMSCORE2e<-0
riskscoredown2<-structure(1:(nrow(geno)), names=rownames(geno))
riskscoreup2<-structure(1:(nrow(geno)), names=rownames(geno))
for (n in 1:(nrow(scores))) {

  rs<-scores$V3[n]
  if (length(which(bim$V2==rs)) == 0) { ## checks lead variant in prs available in test data set, if not skips to next variant
  next
  }

  ## get risk score for lead variant
  beta<-scores[n,"V6"]
  ## check whether A1 allele in test data matches scores A1
  if (bim[which(bim$V2==rs),"V5"]==scores[n,"V4"]) {
    riskscore<-geno[,which(bim$V2==rs)]*beta ## Use beta as is
  } else if (bim[which(bim$V2==rs),"V6"]==scores[n,"V4"]) { ## test A2 in test matches scores A1 
    riskscore<-(geno[,which(bim$V2==rs)]*-1+2)*beta ## invert beta
  } else { ## no match between test and scores. all individuals set to 0
    riskscore[1:(nrow(geno))]<-0
  }
  ## set missing individual genotype calls to have mean riskscore in population
  riskscore[is.na(geno[,which(bim$V2==rs)])]<-mean(riskscore,na.rm=T)

  ## based on rs find closest upstream and downstream variant within ld (0.5 to 0.75) that also exists in test data and sumstats
  pos<-which(ldbim$V2==rs)
  if (length(pos) == 0) { ## checks lead variant available in ld data set, if not adds to running total and moves on to next variant
  ORIGINALSUMSCORE<-ORIGINALSUMSCORE+riskscore
  UPDOGSUMSCORE1a<-UPDOGSUMSCORE1a+riskscore
  UPDOGSUMSCORE1b<-UPDOGSUMSCORE1b+riskscore
  UPDOGSUMSCORE1c<-UPDOGSUMSCORE1c+riskscore
  UPDOGSUMSCORE1d<-UPDOGSUMSCORE1d+riskscore
  UPDOGSUMSCORE1e<-UPDOGSUMSCORE1e+riskscore
  UPDOGSUMSCORE2a<-UPDOGSUMSCORE2a+riskscore
  UPDOGSUMSCORE2b<-UPDOGSUMSCORE2b+riskscore
  UPDOGSUMSCORE2c<-UPDOGSUMSCORE2c+riskscore
  UPDOGSUMSCORE2d<-UPDOGSUMSCORE2d+riskscore
  UPDOGSUMSCORE2e<-UPDOGSUMSCORE2e+riskscore
  next
  }
  start<-min(which(ldbim$V4>as.numeric(ldbim[which(ldbim$V2==rs),4])-window))
  stop<-max(which(ldbim$V4<as.numeric(ldbim[which(ldbim$V2==rs),4])+window))
  suppressWarnings(ldbim$ld[start:stop]<-abs(cor(ld[,start:stop], ld[,pos])))

  ild<-0
  for (i in (pos-1):(start)) { ## loop from 1 variant downstream to beginning of window 
    if (ldbim[i,"ld"] > 0.5 && ldbim[i,"ld"] < 0.75 && (ldbim[i,2] %in% bim$V2) == TRUE && (ldbim[i,2] %in% sumstats$V1) == TRUE) {  ## Find first variant with an ld with lead of >0.5 and <0.75, and is available in the summary stats and test data
      ild<-ldbim[i,"ld"]
      break
    }
  }

  jld<-0
  for (j in (pos+1):(stop)) { ## loop from 1 variant upstream to end of window 
    if (ldbim[j,"ld"] > 0.5 && ldbim[j,"ld"] < 0.75 && (ldbim[j,2] %in% bim$V2) == TRUE && (ldbim[j,2] %in% sumstats$V1) == TRUE && abs(cor(ld[,i],ld[,j])) < 0.9) { ## Find first variant with an ld with lead of >0.5 and <0.75, that is available in the summary stats and test data and isn't in ld (0.9) with the downstream variant
      jld<-ldbim[j,"ld"]
      break
    }
  }

  ## check downstream variant identified
  if (ild > 0) {

    ## Capture causal allele from sumstats scores 
    if (sumstats[which(sumstats$V1==ldbim[i,2]),"V4"] > 0) {
      icausal<-sumstats[which(sumstats$V1==ldbim[i,2]),"V2"]
    } else {
      icausal<-sumstats[which(sumstats$V1==ldbim[i,2]),"V3"]
    }

    ## calculate downstream riskscore
    if (bim[which(bim$V2==ldbim[i,2]),"V5"]==icausal) {  ## If A1 in bim file is causal
      riskscoredown<-geno[,which(bim$V2==ldbim[i,2])]*abs(beta)*ild
      if (mean(riskscore,na.rm=T) > 0) {  ## if lead is causal
        riskscoredown2<-(geno[,which(bim$V2==rs)]-geno[,which(bim$V2==ldbim[i,2])])*-abs(beta)*ild  ##
      } else {   ## if lead is protective
        riskscoredown2<-(geno[,which(bim$V2==rs)]-(geno[,which(bim$V2==ldbim[i,2])]*-1+2))*abs(beta)*ild  ##
      }
    } else if (bim[which(bim$V2==ldbim[i,2]),"V6"]==icausal) { ## If A2 in bim file is causal
      riskscoredown<-geno[,which(bim$V2==ldbim[i,2])]*-abs(beta)*ild
      if (mean(riskscore,na.rm=T) > 0) {  ## if lead is causal
        riskscoredown2<-(geno[,which(bim$V2==rs)]-(geno[,which(bim$V2==ldbim[i,2])]*-1+2))*-abs(beta)*ild  ##
      } else {   ## if lead is protective
        riskscoredown2<-(geno[,which(bim$V2==rs)]-geno[,which(bim$V2==ldbim[i,2])])*abs(beta)*ild  ##
      }
    } else {  ## where test data doesn't match causal allele set riskscoredown to 0
      riskscoredown[1:(nrow(geno))]<-0
      riskscoredown2[1:(nrow(geno))]<-0
    }
    ## set missing individual genotype calls to have riskscore of 0, i.e no correction
    riskscoredown[is.na(geno[,which(bim$V2==ldbim[i,2])])]<-0
    riskscoredown2[is.na(geno[,which(bim$V2==ldbim[i,2])])]<-0
  } else { ## if no downstream variant found set riskscore downstream to 0
    riskscoredown[1:(nrow(geno))]<-0
    riskscoredown2[1:(nrow(geno))]<-0
  }

  if (jld > 0) {

    ## Capture causal allele from sumstats scores
    if (sumstats[which(sumstats$V1==ldbim[j,2]),"V4"] > 0) {
      jcausal<-sumstats[which(sumstats$V1==ldbim[j,2]),"V2"]
    } else {
      jcausal<-sumstats[which(sumstats$V1==ldbim[j,2]),"V3"]
    }
    ## calculate upstream riskscore
    if (bim[which(bim$V2==ldbim[j,2]),"V5"]==jcausal) {
      riskscoreup<-geno[,which(bim$V2==ldbim[j,2])]*abs(beta)*jld
      if (mean(riskscore,na.rm=T) > 0) {  ## if lead is causal
        riskscoreup2<-(geno[,which(bim$V2==rs)]-geno[,which(bim$V2==ldbim[j,2])])*-abs(beta)*jld  ##
      } else {   ## if lead is protective
        riskscoreup2<-(geno[,which(bim$V2==rs)]-(geno[,which(bim$V2==ldbim[j,2])]*-1+2))*abs(beta)*jld  ##
      }
    } else if (bim[which(bim$V2==ldbim[j,2]),"V6"]==jcausal) {
      riskscoreup<-geno[,which(bim$V2==ldbim[j,2])]*-abs(beta)*jld
      if (mean(riskscore,na.rm=T) > 0) {  ## if lead is causal
        riskscoreup2<-(geno[,which(bim$V2==rs)]-(geno[,which(bim$V2==ldbim[j,2])]*-1+2))*-abs(beta)*jld  ##
      } else {   ## if lead is protective
        riskscoreup2<-(geno[,which(bim$V2==rs)]-geno[,which(bim$V2==ldbim[j,2])])*abs(beta)*jld  ##
      }
    } else { ##  where test data doesn't match causal allele set riskscoredown to 0
      riskscoreup[1:(nrow(geno))]<-0
      riskscoreup2[1:(nrow(geno))]<-0
    }
    ## set missing genotype calls to have riskscore of 0
    riskscoreup[is.na(geno[,which(bim$V2==ldbim[j,2])])]<-0
    riskscoreup2[is.na(geno[,which(bim$V2==ldbim[j,2])])]<-0
  } else { ## if no upstream variant found set riskscore upstream to 0
    riskscoreup[1:(nrow(geno))]<-0
    riskscoreup2[1:(nrow(geno))]<-0
  }

  ## set riskscoredown2 and riskscoreup2 to 0 when lead variant is missing
  riskscoredown2[is.na(geno[,which(bim$V2==rs)])]<-0
  riskscoreup2[is.na(geno[,which(bim$V2==rs)])]<-0

  ## running tally of original risk score and updog riskscore
  ORIGINALSUMSCORE<-ORIGINALSUMSCORE+riskscore
  UPDOGSUMSCORE1a<-UPDOGSUMSCORE1a+riskscore+0.025*(riskscoredown+riskscoreup)
  UPDOGSUMSCORE1b<-UPDOGSUMSCORE1b+riskscore+0.05*(riskscoredown+riskscoreup)
  UPDOGSUMSCORE1c<-UPDOGSUMSCORE1c+riskscore+0.125*(riskscoredown+riskscoreup)
  UPDOGSUMSCORE1d<-UPDOGSUMSCORE1d+riskscore+0.25*(riskscoredown+riskscoreup)
  UPDOGSUMSCORE1e<-UPDOGSUMSCORE1e+riskscore+0.5*(riskscoredown+riskscoreup)
  UPDOGSUMSCORE2a<-UPDOGSUMSCORE2a+riskscore+0.025*(riskscoredown2+riskscoreup2)
  UPDOGSUMSCORE2b<-UPDOGSUMSCORE2b+riskscore+0.05*(riskscoredown2+riskscoreup2)
  UPDOGSUMSCORE2c<-UPDOGSUMSCORE2c+riskscore+0.125*(riskscoredown2+riskscoreup2)
  UPDOGSUMSCORE2d<-UPDOGSUMSCORE2d+riskscore+0.25*(riskscoredown2+riskscoreup2)
  UPDOGSUMSCORE2e<-UPDOGSUMSCORE2e+riskscore+0.5*(riskscoredown2+riskscoreup2)

}

output<-cbind(fam[,c(1,2,6)],ORIGINALSUMSCORE,UPDOGSUMSCORE1a,UPDOGSUMSCORE1b,UPDOGSUMSCORE1c,UPDOGSUMSCORE1d,UPDOGSUMSCORE1e,UPDOGSUMSCORE2a,UPDOGSUMSCORE2b,UPDOGSUMSCORE2c,UPDOGSUMSCORE2d,UPDOGSUMSCORE2e)
colnames(output)[1:3]<-c("FID","ID","PHENO")

## write out risk score
write.table(output,paste0("temp/scoreoutput_",name,"_chr",chr,"_",sprintf("%03d",chunk),".txt"),col.names=T,row.names=F,quote=FALSE)