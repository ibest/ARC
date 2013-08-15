library(ShortRead)

for(f in dir(pattern="finished_*")){
  cat(f)
  cat("\n")
  if(file.exists(paste(f,"/contigs.fasta", sep=""))){
      fna = readDNAStringSet(paste(f,"/contigs.fasta", sep=""), "fasta")
      if(length(fna) > 0){
        t = table(table(sapply(strsplit(names(fna), "_:_", fixed=T), '[', 2)))
        cat(names(t))
        cat("\n")
        cat(t)
        cat("\n")
        }else{
      cat(0)
      }
    }else
      {
        cat(0)
      }
    cat('\n')

}


