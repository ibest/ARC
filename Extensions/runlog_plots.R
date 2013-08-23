### v1  --  Split running time ####
m = function(t){
    t = as.numeric(t)
    t[1]*60*60 + t[2]*60 + t[3]
}

#while(TRUE){

## Splits
system("grep 'Split' log.txt > splits")
splits = read.table("splits", fill=T, as.is=T)
splits = splits[splits$V16 == "seconds" & splits$V5 == "Sample:",]

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
assemblies = assemblies[assemblies$V15 == "seconds" & assemblies$V5 == "Sample:",]

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
finished = read.table("finished", fill=T)
finished = finished[finished$V20  == "done" & finished$V5 == "Sample", ]

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
killed = read.table("killed", fill=T)
killed = killed[killed$V5  == "Sample:" & killed$V15 == "seconds", ]

#Figure out finished per second:
killed_p_s = table( paste(killed$V1, sapply(strsplit(killed$V2, ","), '[', 1), sep="_") )
ts = sapply(strsplit(sapply(strsplit(names(killed_p_s), "_"), '[', 2), ":" ),m)
days = as.numeric(sapply(strsplit(sapply(strsplit(names(killed_p_s), "_"), '[', 1), "-" ), '[', 3))
days = days - mdays
killed_p_s_times = ts + days*24*60*60
kiled_p_s_times  = killed_p_s_times - zero

## repetitive:
system('grep -E  "repetitive" log.txt > repeat')
repeats = read.table("repeat", fill=T)
repeats = repeats[repeats$V5  == "Sample" & repeats$V18 == "done", ]

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

par(mfrow=c(1,6))
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
text(x=1, y=310000, labels=paste("Total killed:", sum(killed_p_s), sep=" "), adj=c(0,0))
text(x=1, y=300000, labels=paste("Total completed:", sum(finished_p_s), sep=" "), adj=c(0,0))
text(x=1, y=290000, labels=paste("Total repeats:", sum(repeats_p_s), sep=" "), adj=c(0,0))


dim(splits)[1] - dim(assemblies)[1] - dim(killed)[1]

Sys.sleep(20)
}




########################### ------------------- ##################
# Check splits vs finished
sp_uniq = paste(splits$V6, splits$V8, splits$V10, sep="_")
as_uniq = paste(assemblies$V6, assemblies$V8, assemblies$V10, sep="_")
ki_uniq = paste(killed$V6, killed$V8, killed$V10, sep="_")

(sp_uniq %in% as_uniq | sp_uniq %in% ki_uniq)
sp_uniq[!(sp_uniq %in% as_uniq | sp_uniq %in% ki_uniq))]

















### v2  --  Split running time ####
m = function(t){
    t = as.numeric(t)
    t[1]*60*60 + t[2]*60 + t[3]
}

while(TRUE){

system("grep 'Split' log.txt > splits")
splits = read.table("splits", fill=T)
#splits = na.omit(dat)
splits = splits[splits$V12 == "Split",]
splits_p_s = table( sapply(strsplit(splits$V2, ","), '[', 1) )

splits_p_s_times = sapply(strsplit(names(splits_p_s), ':'), m)
zero = splits_p_s_times[1]
splits_p_s_times = splits_p_s_times - zero

system("grep 'Assembly finished' log.txt > assemblies")
assemblies = read.table("assemblies", fill=T)
#assemblies = na.omit(assemblies)
assemblies = assemblies[assemblies$V12 == "Assembly",]
assemblies_p_s = table( sapply(strsplit(assemblies$V2, ","), '[', 1) )
assemblies_p_s_times = sapply(strsplit(names(assemblies_p_s), ':'), m)
assemblies_p_s_times = assemblies_p_s_times - zero

par(mfrow=c(1,5))
plot(splits$V16, main="Splitting times", xlab="split", ylab="seconds", col=as.factor(splits$V10), pch=".")
plot(splits_p_s, main="Splits per second", xlab="time", ylab="Splits/Second", type='p', pch=".")

plot(assemblies$V15, main="Assembly times", xlab="assembly", ylab="seconds", col=as.factor(splits$V10), pch=".")
plot(assemblies_p_s, main="Assemblies per second", xlab="time", ylab="Assemblies/Second", type='p', pch=".")

plot(x=splits_p_s_times, y=cumsum(splits_p_s), pch=".", col="red", cex=2, xlab="Time (s)", ylab = "Operations performed",
    xlim=range(c(splits_p_s_times, assemblies_p_s_times)), ylim=c(0, sum(splits_p_s)))
points(x=assemblies_p_s_times, y=cumsum(assemblies_p_s), pch=".", col="black", cex=2)
legend("topleft", col=c("red","black"), pch=20, legend=c("Splits", "Assemblies"))


Sys.sleep(20)
}
















### Assemblies per second ####
system("grep Split log.txt > tmp2")
splits = read.table("tmp", fill=T)
splits = na.omit(splits)

splits_p_s = table( sapply(strsplit(splits$V2, ","), '[', 1) )
splits_p_m = table( substr(splits$V2, 1, 5))
}


### Assembly running time ####
system("grep 'Assembly finished' log.txt > tmp")
dat = read.table("tmp", fill=T)
dat = na.omit(dat)
#plot(dat$V14, pch=".", cex=2, main="Newbler assembly time",
#    xlab="Consecuitive Assembly", ylab="Time (s)", col=as.factor(dat$V8))

plot(dat$V14, pch=".", cex=2, main="Newbler assembly time",
    xlab="Consecuitive Assembly", ylab="Time (s)", col=dat$V10)

#legend("topright", col=unique(as.factor(dat$V8)), pch=20, legend=unique(dat$V8), cex=.4)
legend(title="Iteration", "topright", col=unique(as.factor(dat$V10)), pch=20, legend=unique(dat$V10), cex=.8)

m = 1:dim(dat)[1]
l = lm(dat$V14~m)
abline(l, col="black", lwd=2)

w = median(dat$V14)/dat$V14
w[dat$V14 < 1] = 1
ss = smooth.spline(m,  w=w, dat$V14, spar=0.35)
lines(ss, col="orange", lwd=2)

boxplot(dat$V14~dat$V10)

#lo = loess(m~dat$V12)
#lines(predict(lo), col='red', lwd=2)
#lo = loess(dat$V12~m)
#lines(predict(lo), col='red', lwd=2)

# Get contigs lengths:

lens = fasta.info("100_targets_genes.txt")
targets = sapply(strsplit(names(lens), "_:_", fixed=T), '[', 2)
tlenghts = tapply(lens, INDEX=targets, sum)
dat2 = cbind(dat, tlengths=tlenghts[match(dat$V8, names(tlenghts))])
boxplot(dat$V12~dat$tlengths)
cor(dat2$V12,dat2$tlengths)



### Failed assemblies ####
system("grep -i kill spades_log.txt > tmp")
dat = read.table("tmp")

barplot(table(dat$V8))

plot(dat$V12, pch=".", cex=2, main="Spades assembly time",
    xlab="Consecuitive Assembly", ylab="Time (s)", col=as.factor(dat$V8))
#legend("topright", col=unique(as.factor(dat$V8)), pch=20, legend=unique(dat$V8), cex=.4)

m = 1:dim(dat)[1]
l = lm(dat$V12~m+m^2)
abline(l, col="black", lwd=2)

ss = smooth.spline(m, dat$V12, spar=0.35)
lines(ss, col="red", lwd=2)
