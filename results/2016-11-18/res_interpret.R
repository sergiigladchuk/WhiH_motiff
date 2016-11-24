
# all peak analysis
fimo <- read.csv('2016-11-18/all_peaks/fimo.txt', sep='\t')
#View(fimo)


fimo$logs <- -log(fimo$p.value)
# add matches from MCAST

mcast <- read.csv('2016-11-18/all_peaks/mcast.txt', sep='\t')
#View(mcast)


fimo$mcast <- match(fimo$start,mcast$start)
fimo$mcast[!is.na(fimo$mcast)] <- fimo$logs[!is.na(fimo$mcast)]

# add matches from Centrimo



#plot all results


plot(fimo$logs, type = 'l', main="-log of p-value vs. index", 
     ylab = "-log of p-value", xlab = "# of matches")
points(fimo$mcast, col='red')
