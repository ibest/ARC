#!/usr/bin/env Rscript

cat("\nGenerate memory usage plots from data collected by profilemem.R.\n")
cat("\nBy default plots are drawn to a PDF, supply a file-name argument to create a file instead.\n")
cat("\nAlternatively you can source this file from within R to get plots to an x11() device\n")

cat("Usage: plot_memprofile.R <filename>\n\n")

plottofile = FALSE

if(length(commandArgs(trailingOnly=T)) == 1){
  outfname = commandArgs(trailingOnly=T)[1]
  plottofile = TRUE
  print(paste("Plotting to", outfname))
}

if(plottofile){
  png(outfname, width=10, height=5, units='in')
  }else{
    dev.new(width=10, height=5)
}

time_data = read.table("time_data.tsv", fill=T, as.is=T, header=T)
time_data$Mem = as.numeric(time_data$Mem)
time_data$CPU = as.numeric(time_data$CPU)
time_data = na.omit(time_data)

#t = table(time_data$PID)
#idx = names(t[t>1])
#time_data2 = time_data[time_data$PID %in% idx,]
#time_data2$time = time_data2$time - min(time_data2$time)

#png(file="memory_usage_over_time.png", width=2000, height=1000)
  layout(mat = matrix(c(1,2,3,4), nrow=2, byrow=T), heights=c(3,5), widths=1)
  cols = rainbow(length(unique(time_data$PID)))
  names(cols) = unique(time_data$PID)

  par(mar=c(0,1,1,1))
  plot(NA, xlim=range(0,25), ylim=range(0,25), xaxt='n', yaxt='n', ylab="", xlab="")
  #legend("top", col=cols, legend=unique(time_data2$PID), pch=20, ncol=30, title="ProcessID", bty='n', cex=.6)

  par(mar=c(2,4,1,1))
  total_mem = tapply(as.numeric(time_data$Mem), INDEX=time_data$time, FUN=sum)
  total_mem_time = as.numeric(names(total_mem))
  plot(y=total_mem/1024/1024, x=total_mem_time, type='l', ylab="Total memory (Gb)")


  par(mar=c(4,4,2,3))
  plot(NA, xlim=range(time_data$time), ylim=range(time_data$Mem/1024), ylab="Memory (MB)", xlab="Seconds")
  #for(p in unique(time_data2$PID)){
  #  col = cols[as.character(p)]
  #  lines((Mem/1024)~time, time_data2[time_data2$PID == p,], col=col, type="o", pch=".")
  #}
  points(y=time_data$Mem/1024, x=time_data$time, col=cols[as.character(time_data$PID)], pch=".", cex=2)


  par(mar=c(4,4,2,3))
  plot(NA, xlim=range(time_data$time), ylim=range(time_data$CPU), ylab="Percent CPU", xlab="Seconds")
  #for(p in unique(time_data2$PID)){
  #  col = cols[as.character(p)]
  #  lines((CPU)~time, time_data2[time_data2$PID == p,], col=col, type="o", pch=".")
  #}
  points(y=time_data$CPU, x=time_data$time, col=cols[as.character(time_data$PID)], pch=".", cex=2)

cat("Plots created...\n")

if(plottofile){
  dev.off()
}
