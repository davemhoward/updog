chunkfiles <- list.files(pattern="chunkscores_chr*")
#length(chunkfiles)
scores<-paste0(gsub("chunkscores", "scoreoutput", chunkfiles),".txt")

p<-0
for (n in 1:(length(chunkfiles))) {
  if ((file.exists(scores[n])) == FALSE) {
    p<-1
    print(paste0(scores[n]," is missing"))
  }
}

if (p == 1) {
  stop("need to rerun missing file(s). try extending runtime by entering the following on the command line for each missing file replacing the square brackets:\n qsub -t [chunk]-[chunk] -l h_rt=4:00:00 ./updog.sh [chr]")
}

set1<-read.table(scores[1],header=T,sep="")

for (i in 2:(length(scores))){
#for (i in 2:2) {
  set2<-read.table(scores[i],header=T,sep="")
  set1[,4:(ncol(set1))]<-set1[,4:(ncol(set1))] + set2[,4:(ncol(set2))]
}

write.table(set1,"scoreoutput_genomewide.txt",quote=F,col.names=T,row.names=F)
