
# all peak analysis
fimo <- read.csv('2016-12-05/all/fimo.txt', sep='\t')
#View(fimo)


fimo$logs <- -log10(fimo$p.value)
# add matches from MCAST

# add ranges from all chip
chip_all <- read.csv('2016-11-18/all_peaks/centrimo.txt')

# add ranges from top chip
chip_top <- read.csv

fimo$centrimo_full_match <- c(NA)
totMathes <- 0
rangeExtantion <- 100
for (count1 in 1:length(centrimo$Centrimo_Reg_Matches)) {
  for (count2 in 1:length(fimo$logs)) {
    if (centrimo$Centrimo_Reg_Matches[count1] >= (fimo$start[count2] - rangeExtantion)
        && centrimo$Centrimo_Reg_Matches[count1] <= (fimo$stop[count2] + rangeExtantion)){
      fimo$centrimo_full_match[count2] = fimo$logs[count2]
      
      totMathes = totMathes + 1
      break;
    }
  }
  
}
cat(totMathes, "CentriMo regions in all FIMO out of", length(centrimo$Centrimo_Reg_Matches), 
    "with", rangeExtantion ,"range extention value")



#plot all results


plot(fimo$logs, type = 'l', main="-log of p-value vs. index", 
     ylab = "-log of p-value", xlab = "index in FIMO")
points(fimo$mcast, col='red')
points(fimo$centrimo_full_match, col='blue',pch=3)
legend(x='topright',legend = c('FIMO','MCAST','CentriMo'), lty=c(1,0,0),pch = c(NA,1,3), col=c('black','red','blue'))
