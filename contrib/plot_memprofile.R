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

total_CPU = tapply(as.numeric(time_data$CPU), INDEX=time_data$time, FUN=sum)
total_CPU_time = as.numeric(names(total_CPU))
total_mem = tapply(as.numeric(time_data$Mem), INDEX=time_data$time, FUN=sum)
total_mem_time = as.numeric(names(total_mem))

options(scipen=2)
layout(mat = matrix(c(1,2,3,4), nrow=2, byrow=T), heights=c(3,5), widths=1)
cols = rainbow(length(unique(time_data$PID)))
names(cols) = unique(time_data$PID)

par(mar=c(2,5,1,1))
plot(y=total_CPU, x=total_CPU_time, type='l', ylab="Total % CPU", cex.axis=1.2, cex.lab=1.2)
text(y=max(total_CPU)/5.5, x=max(time_data$time)/2,
     paste("Running time", signif((max(time_data$time) - min(time_data$time))/60/60, 4), "hrs\n",
     "Max cpu:", signif(max(total_CPU), 4), "%\n",
     "Average cpu:", signif(mean(total_CPU), 4), "%" ))


par(mar=c(2,4,1,1))
plot(y=total_mem/1024/1024, x=total_mem_time, type='l', ylab="Total memory (Gb)", cex.axis=1.2, cex.lab=1.2)
text(y=max(total_mem)/1024/1024/2.5, x=max(time_data$time)/2,
     paste("Max memory", signif(max(total_mem)/1024/1024, 4), "Gb\n",
     "Average memory:", signif(mean(total_mem)/1024/1024, 4), "Gb\n"))

par(mar=c(4,5,2,1))
plot(NA, xlim=range(time_data$time), ylim=range(time_data$Mem/1024), ylab="Per-Process Memory (MB)", xlab="Seconds", cex.axis=1.2, cex.lab=1.2)
points(y=time_data$Mem/1024, x=time_data$time, col=cols[as.character(time_data$PID)], pch=".", cex=2)

par(mar=c(4,4,2,1))
plot(NA, xlim=range(time_data$time), ylim=range(time_data$CPU), ylab="Per-Process Percent CPU", xlab="Seconds", cex.axis=1.2, cex.lab=1.2)
points(y=time_data$CPU, x=time_data$time, col=cols[as.character(time_data$PID)], pch=".", cex=2)

cat("Plots created...\n")

if(plottofile){
  dev.off()
}
