### v1  --  Split running time ####
m = function(t){
    t = as.numeric(t)
    t[1]*60*60 + t[2]*60 + t[3]
}

## Splits
system("grep 'Split' log.txt > splits")
splits = read.table("splits", fill=T, as.is=T)
splits = splits[splits$V16 == "seconds" & splits$V5 == "Sample:" & substr(splits$V1, 1,1) == '[',]

#Figure out split times in seconds:
splits_times = sapply(strsplit( sapply(strsplit(splits$V2, ","), '[', 1), ':'), m)
day = as.numeric(sapply(strsplit(splits$V1, "-"), '[', 3))
day = day - min(day)
splits_times = splits_times + day*24*60*60
zero = min(splits_times)
splits_times = splits_times - zero

#Figure out splits per second:
splits_p_s = table( paste(splits$V1, sapply(strsplit(splits$V2, ","), '[', 1), sep="_") )
ts = sapply(strsplit(sapply(strsplit(names(splits_p_s), "_"), '[', 2), ":" ),m)
days = as.numeric(sapply(strsplit(sapply(strsplit(names(splits_p_s), "_"), '[', 1), "-" ), '[', 3))
days = days - min(days)
splits_p_s_times = ts + days*24*60*60
splits_p_s_times  = splits_p_s_times - zero


## Assemblies
system("grep 'Assembly finished' log.txt > assemblies")
assemblies = read.table("assemblies", fill=T, as.is=T)
assemblies = assemblies[assemblies$V15 == "seconds" & assemblies$V5 == "Sample:" & substr(assemblies$V1, 1,1) == '[',]

#Figure out assemblies times in seconds:
ts = sapply(strsplit( sapply(strsplit(assemblies$V2, ","), '[', 1), ':'), m)
days = as.numeric(sapply(strsplit(assemblies$V1, "-"), '[', 3))
days = days - min(days)
assembly_times = ts + days*24*60*60
assembly_times = assembly_times - zero

#Figure out assemblies per second
assemblies_p_s = table( paste(assemblies$V1, sapply(strsplit(assemblies$V2, ","), '[', 1), sep="_") )
ts = sapply(strsplit(sapply(strsplit(names(assemblies_p_s), "_"), '[', 2), ":" ),m)
days = as.numeric(sapply(strsplit(sapply(strsplit(names(assemblies_p_s), "_"), '[', 1), "-" ), '[', 3))
mdays = min(days)
days = days - mdays
assemblies_p_s_times = ts + days*24*60*60
assemblies_p_s_times  = assemblies_p_s_times - zero

## Finished:
system('grep "did not incorporate any more reads" log.txt > finished')
finished = read.table("finished", fill=T, as.is=T)
finished = finished[finished$V20  == "done" & finished$V5 == "Sample" & substr(finished$V1, 1,1) == '[', ]

#Figure out finished per second:
ts = sapply(strsplit( sapply(strsplit(finished$V2, ","), '[', 1), ':'), m)
days = as.numeric(sapply(strsplit(finished$V1, "-"), '[', 3))
days = days - mdays
finished_times = ts + days*24*60*60
finished_times = assembly_times - zero

#Figure out finished per second:
finished_p_s = table( paste(finished$V1, sapply(strsplit(finished$V2, ","), '[', 1), sep="_") )
ts = sapply(strsplit(sapply(strsplit(names(finished_p_s), "_"), '[', 2), ":" ),m)
days = as.numeric(sapply(strsplit(sapply(strsplit(names(finished_p_s), "_"), '[', 1), "-" ), '[', 3))
days = days - mdays
finished_p_s_times = ts + days*24*60*60
finished_p_s_times  = finished_p_s_times - zero

## killed:
system('grep "Assembly killed" log.txt > killed')
killed = read.table("killed", fill=T, as.is=T)
killed = killed[killed$V5  == "Sample:" & killed$V15 == "seconds" & substr(killed$V1, 1,1) == '[', ]

#Figure out finished per second:
killed_p_s = table( paste(killed$V1, sapply(strsplit(killed$V2, ","), '[', 1), sep="_") )
ts = sapply(strsplit(sapply(strsplit(names(killed_p_s), "_"), '[', 2), ":" ),m)
days = as.numeric(sapply(strsplit(sapply(strsplit(names(killed_p_s), "_"), '[', 1), "-" ), '[', 3))
days = days - mdays
killed_p_s_times = ts + days*24*60*60
kiled_p_s_times  = killed_p_s_times - zero

## repetitive:
system('grep -E  "repetitive" log.txt > repeat')
repeats = read.table("repeat", fill=T, as.is=T)
repeats = repeats[repeats$V5  == "Sample" & repeats$V18 == "done" & substr(repeats$V1, 1,1) == '[', ]

#Figure out finished per second:
repeats_p_s = table( paste(repeats$V1, sapply(strsplit(repeats$V2, ","), '[', 1), sep="_") )
ts = sapply(strsplit(sapply(strsplit(names(repeats_p_s), "_"), '[', 2), ":" ),m)
days = as.numeric(sapply(strsplit(sapply(strsplit(names(repeats_p_s), "_"), '[', 1), "-" ), '[', 3))
days = days - mdays
repeats_p_s_times = ts + days*24*60*60
repeats_p_s_times  = repeats_p_s_times - zero



## Make plots
xlim = range(splits_times)

d = 3600
xlab = "time (hr)"

#par(mfrow=c(2,3))
layout(mat=matrix(c(1,2,3,4,5,6), nrow=2, byrow=F), heights=c(1,1))
plot(splits_times/d, as.numeric(splits$V15), main="Splitting times", xlab=xlab, ylab="seconds", col=as.factor(splits$V10), pch=".", xlim=xlim/d)
plot(y=as.numeric(splits_p_s), x=splits_p_s_times/d, main="Splits per second", xlab=xlab, ylab="Splits/Second", type='p', pch=".", xlim=xlim/d)

plot(assembly_times/d, as.numeric(assemblies$V14), main="Assembly times", xlab=xlab, ylab="seconds", col=as.factor(assemblies$V10), pch=".", xlim=xlim/d)
plot(assemblies_p_s_times/d, y=as.numeric(assemblies_p_s), main="Assemblies per second", xlab=xlab, ylab="Assemblies/Second", type='p', pch=".", xlim=xlim/d)

plot(x=splits_p_s_times/d, y=cumsum(splits_p_s), pch=".", col="red", cex=2, xlab=xlab, ylab = "Operations performed",
    xlim=xlim/d, ylim=c(0, sum(splits_p_s)), main="Total Splits and Assemblies")
points(x=assemblies_p_s_times/d, y=cumsum(assemblies_p_s), pch=".", col="black", cex=2)
legend("topleft", col=c("red","black"), pch=20, legend=c("Splits", "Assemblies"))

plot(x=finished_p_s_times/d, y=cumsum(finished_p_s), pch=".", col="black", cex=2, xlab=xlab, ylab = "Total complete",
    xlim=xlim/d, ylim=c(0, sum(finished_p_s)), main="Completed Targets")
points(x=killed_p_s_times/d, y=cumsum(killed_p_s), pch=".", col="red", cex=2)
points(x=repeats_p_s_times/d, y=cumsum(repeats_p_s), pch=".", col="green", cex=2)
legend("topleft", col=c("black","red", "green"), pch=20, legend=c("Completed", "Killed", "Repeatitive"))
text(x=max(xlim/d) * .01, y=sum(finished_p_s) * .75, labels=paste("Total completed:", sum(finished_p_s), sep=" "), adj=c(0,0))
text(x=max(xlim/d) * .01, y=sum(finished_p_s) * .7, labels=paste("Total killed:", sum(killed_p_s), sep=" "), adj=c(0,0))
text(x=max(xlim/d) * .01, y=sum(finished_p_s) * .65, labels=paste("Total repeats detected:", sum(repeats_p_s), sep=" "), adj=c(0,0))




