args <- commandArgs(trailingOnly = TRUE)

name <- args[1]

setwd("temp")

chunkname<-paste0("chunkscores_",name,"_chr*")
chunkfiles <- list.files(pattern=chunkname)
scores<-paste0(gsub("chunkscores", "scoreoutput", chunkfiles),".txt")

## sapply(strsplit(sapply(strsplit(scores[n], "_chr"), "[[" , 2), "_"), "[[" , 1) ## this captures the missing chromosome from scores
## as.numeric(substr(scores[n],nchar(scores[n])-6,nchar(scores[n])-4)) ## this captures the missing chunk from scores

p<-0
for (n in 1:(length(chunkfiles))) {
  if ((file.exists(scores[n])) == FALSE) {
    p<-p+1
    cat(paste0(  scores[n]," is missing \n"))
    if (p == 1) {
      sink("../resubmitjobs",append=FALSE)
    } else {
      sink("../resubmitjobs",append=TRUE)
    }
    cat("qsub -t 1- -l h_rt=4:00:00 -o logs/ -e logs/ -l h_vmem=16G -N updog -cwd ./updog.sh args go here \n")
    sink()
  }
}

if (p >= 1) {
  ## add additional job to the resubmitjobs file that reruns this merge at the end
  stop("To rerun missing chunks type ./resubmitjobs on the command line") ## COULD CREATE A -z flag for testlaunch to just rerun the merge afterwards
}

set1<-read.table(scores[1],header=T,sep="")

for (i in 2:(length(scores))){
  set2<-read.table(scores[i],header=T,sep="")
  set1[,4:(ncol(set1))]<-set1[,4:(ncol(set1))] + set2[,4:(ncol(set2))]
}

setwd("..")   
ADD NAME TO OUTPUT
write.table(set1,"scoreoutput_genomewide.txt",quote=F,col.names=T,row.names=F)
