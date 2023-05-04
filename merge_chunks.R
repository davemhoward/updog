# supplied arguments: $1=$testloc, $2=$testtype, $3=$ldloc, $4=$ldtype, $5=$sumstats, $6=$scores
# $7=$plinkloc, $8=$rloc, $9=$outname, $10=sumstatsOR, $11=jobargs, $12=weighting

args <- commandArgs(trailingOnly = TRUE)
name <- args[9]

outargs<-paste(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],args[10],args[11],args[12])

setwd("temp")
chunkname<-paste0("chunkscores_",name,"_chr*")
chunkfiles <- list.files(pattern=chunkname)
scores<-paste0(gsub("chunkscores", "scoreoutput", chunkfiles),".txt")

p<-0
for (n in 1:(length(chunkfiles))) {
  if ((file.exists(scores[n])) == FALSE) {
    p<-p+1
    cat(paste0(  scores[n]," is missing \n"))
    chr<-sapply(strsplit(sapply(strsplit(scores[n], "_chr"), "[[" , 2), "_"), "[[" , 1)
    chunk<-as.numeric(substr(scores[n],nchar(scores[n])-6,nchar(scores[n])-4))

    if (p == 1) {
      sink(paste0("../resubmitjobs_",name),append=FALSE)
    } else {
      sink(paste0("../resubmitjobs_",name),append=TRUE)
    }
    sbatchargs<-paste0("--export=ALL,i=",chr,",testloc=",args[1],",testtype=",args[2],",ldloc=",args[3],",ldtype=",args[4],",sumstats=",args[5],",scores=",args[6],",plinkloc=",args[7],",rloc=",args[8],",outname=",args[9],",sumstatsOR=",args[10],",jobargs=",args[11],",weighting=",args[12])
    cat(paste0("sbatch --array=",chunk,"-",chunk," --time=0-8:00:00 --output=logs/resubmit.chr",chr,".",chunk,".txt ",args[11]," --mem=32G --job-name=updog_resubmit ",sbatchargs," ./chunk.sh\n"))
    sink()

  }
}

if (p >= 1) {

  sink(paste0("../resubmitjobs_",name),append=TRUE)
  mergeargs<-paste0("--export=ALL,testloc=",args[1],",testtype=",args[2],",ldloc=",args[3],",ldtype=",args[4],",sumstats=",args[5],",scores=",args[6],",plinkloc=",args[7],",rloc=",args[8],",outname=",args[9],",sumstatsOR=",args[10],",jobargs=",args[11],",weighting=",args[12])
  cat(paste0("sbatch --time=0-8:00:00 --output=logs/resubmit_merge_chunks.txt ",args[11]," --mem=32G --job-name=updog_resubmit --dependency=singleton ",mergeargs," ./merge_chunks.sh\n"))
  sink()

  cat(paste0("  To rerun missing chunk(s) and attempt a remerge type ./resubmitjobs_",name," on the command line depending on your job scheduler\n"))
  stop()
}

set1<-read.table(scores[1],header=T,sep="")

for (i in 2:(length(scores))){
  set2<-read.table(scores[i],header=T,sep="")
  set1[,4:(ncol(set1))]<-set1[,4:(ncol(set1))] + set2[,4:(ncol(set2))]
}

setwd("..")   
write.table(set1,paste0("genomewidescores_",name,".txt"),quote=F,col.names=T,row.names=F)
