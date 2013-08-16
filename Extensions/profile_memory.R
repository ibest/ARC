user = "shunter"

#dat = read.table("memuse.tsv")
#colnames(dat) = c("User","Mem","PID","ARGS")
#cbind(dat, time = rep(proc.time()[3], dim(dat)[1]))

time_data = data.frame()

while(TRUE){
  if(user == ""){
    system('ps xo euser,rss,pid,%cpu,args | grep ARC | grep python | grep -v "grep" > memuse.tsv')
  }else{
    system(paste('ps xo euser,rss,pid,%cpu,args | grep ARC | grep python | grep', user,  '| grep -v "grep" > memuse.tsv'))
  }
  if(file.exists("memuse.tsv")){
    dat = read.table("memuse.tsv")
    colnames(dat) = c("User","Mem","PID","CPU","ARGS")
    dat = cbind(dat, time = rep(proc.time()[3], dim(dat)[1]))
    time_data = rbind(time_data, dat[,c("Mem","PID","CPU","time")])
    write.table(time_data, "time_data.tsv")
    }
  Sys.sleep(10)
}

while(TRUE){
  time_data = read.table("time_data.tsv", fill=T)
  time_data = na.omit(time_data)
  t = table(time_data$PID)
  idx = names(t[t>2])
  time_data = time_data[time_data$PID %in% idx,]
  
  
  #png(file="memory_usage_over_time.png", width=2000, height=1000)
    layout(mat = matrix(c(1,2,3), nrow=3), heights=c(3,5,5), widths=1)
    cols = rainbow(length(unique(time_data$PID)))
    names(cols) = unique(time_data$PID)

    par(mar=c(0,1,1,1))
    plot(NA, xlim=range(0,25), ylim=range(0,25), xaxt='n', yaxt='n', ylab="", xlab="")
    legend("top", col=cols, legend=unique(time_data$PID), pch=20, ncol=35, title="ProcessID", bty='n', cex=.8)
    
    par(mar=c(1,4,2,3))
    plot(NA, xlim=range(time_data$time), ylim=range(time_data$Mem/1024/1024), ylab="Memory (GB)", xaxt="n", xlab="")
    for(p in unique(time_data$PID)){
      col = cols[as.character(p)]
      lines((Mem/1024/1024)~time, time_data[time_data$PID == p,], col=col, type="o", pch=".")
    }
    #points(y=time_data$Mem/1024, x=time_data$time, col=cols[as.character(time_data$PID)], pch=".", cex=2)
    
    plot(NA, xlim=range(time_data$time), ylim=range(time_data$CPU), ylab="Percent CPU", xlab="Seconds")
    for(p in unique(time_data$PID)){
      col = cols[as.character(p)]
      lines((CPU)~time, time_data[time_data$PID == p,], col=col, type="o", pch=".")
    }
    Sys.sleep(20)
}
#dev.off()

  
  assemblies:
19477
Splits:
167717

assemblies:
21778
Splits:
182230


