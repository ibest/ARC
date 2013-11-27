############################################
library(ShortRead)
library(Biostrings)

libraries = list()

for(f in dir(pattern="finished_*")){
  cat(paste(f, "total completed"))
  cat("\n")
  if(file.exists(paste(f,"/contigs.fasta", sep=""))){
      fna = readDNAStringSet(paste(f,"/contigs.fasta", sep=""), "fasta")
      if(length(fna) > 0){
        t = table(table(sapply(strsplit(names(fna), "_:_", fixed=T), '[', 2)))
        cat(names(t))
        cat("\n")
        cat(t)
        cat("\n")
        libraries[[f]] = table(sapply(strsplit(names(fna), "_:_", fixed=T), '[', 2))
        }else{
      cat(0)
      }
    }else
      {
        cat(0)
      }
    cat('\n')
}

# Determine how many have 1 contig in all:
alltargets = unique( sapply( strsplit(names(unlist(libraries)), ".", fixed=T), '[', 2))

contig_counts = matrix(0, nrow=length(alltargets), ncol=length(libraries))
rownames(contig_counts) = alltargets
colnames(contig_counts) = sapply(strsplit(names(libraries), "_"), '[', 2)

for(library in names(libraries)){
  lib = unlist(strsplit(library, "_"))[2]
  idx = match(names(libraries[[library]]), alltargets)
  contig_counts[idx,lib] = libraries[[library]]
}

# idx3 = order(colSums(contig_counts < 3))
# barplot(colSums(contig_counts < 3)[idx3], main="Number of targets with < 3 contigs")
par(las=2)
idx1 = order(colSums(contig_counts == 1))
barplot(colSums(contig_counts == 1)[idx1], main="Number of targets with 1 contig")


barplot(table(rowSums(contig_counts == 0)), main="Number of targets for which no contig was assembled",
        ylab="Targets", xlab="Number of missing samples")

# Determine how many have 1 contig
barplot(table(rowSums(contig_counts == 1)), ylab="Number of targets", xlab="Number of samples",
        main="Number of targets with a single contig across samples")
