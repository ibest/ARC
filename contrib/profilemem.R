#!/usr/bin/env Rscript

if(length(commandArgs(trailingOnly=T)) == 1){
  user = commandArgs(trailingOnly=T)[1]
  print(user)
} else {
  cat("\nA simple program which uses ps to gather memory and CPU usage for all\n")
  cat("of a user's processes, then write them to a time_data.tsv file.\n")
  cat("The program continues to run until terminated with cntrl+c.\n")
  cat("Usage: profilemem.R <userid>\n\n")
  q()
}


#dat = read.table("memuse.tsv")
#colnames(dat) = c("User","Mem","PID","ARGS")
#cbind(dat, time = rep(proc.time()[3], dim(dat)[1]))

time_data = data.frame()
if(file.exists("time_data.tsv")){
  cat("Error, time_data.tsv already exists, please move or delete it and re-run this program.")
  q()
}

#file.remove("time_data.tsv")
cat("Mem\tPID\tCPU\ttime\n", file="time_data.tsv")

while(TRUE){
  if(user == ""){
    #system('ps xo euser,rss,pid,%cpu,args | grep ARC | grep python | grep -v "grep" > memuse.tsv')
    system('ps xo euser,rss,pid,%cpu | grep -Ev "grep|ps" > memuse.tsv')
  }else{
    #system(paste('ps xo euser,rss,pid,%cpu,args | grep ARC | grep python | grep', user,  '| grep -v "grep" > memuse.tsv'))
    system(paste('ps xo euser,rss,pid,%cpu | grep', user,  '| grep -Ev "grep|ps" > memuse.tsv'))
  }
  if(file.exists("memuse.tsv")){
    dat = read.table("memuse.tsv", fill=T, header=F, row.names=NULL, quote="", comment.char="")
    #colnames(dat) = c("User","Mem","PID","CPU","ARGS")
    colnames(dat) = c("User","Mem","PID","CPU")
    dat = cbind(dat, time = rep(proc.time()[3], dim(dat)[1]))
    #time_data = rbind(time_data, dat[,c("Mem","PID","CPU","time")])
    time_data = dat[,c("Mem","PID","CPU","time")]
    write.table(time_data, "time_data.tsv", append=T, col.names=F, row.names=F)
    }
  Sys.sleep(10)
}

