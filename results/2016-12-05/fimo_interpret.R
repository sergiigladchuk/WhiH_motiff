
# all peak analysis
fimo <- read.csv('2016-12-15/memechip_10000rand_0bg_pal/fimo.txt', sep='\t')
#View(fimo)

# palindrom correction
 fimo <- fimo[-(seq(2,to=nrow(fimo),by=2)),]


fimo$logs <- -log10(fimo$p.value)
# add matches from MCAST

# add ranges from all chip
chip_all <- read.csv('2016-12-05/all_ranges.txt', header = FALSE)

# add ranges from top chip
chip_top <- read.csv('2016-12-05/top_ranges.txt', header = FALSE)

fimo$all_match <- c(NA)
fimo$top_match <- c(NA)
allMatches <- 0

rangeExtantion <- 100 # half and +25 of 100, which was the range of MEME re-discovery
for (count1 in 1:length(chip_all$V1)) {
  for (count2 in 1:length(fimo$logs)) {
    if (chip_all$V1[count1] >= (fimo$start[count2] - rangeExtantion)
        && chip_all$V1[count1] <= (fimo$stop[count2] + rangeExtantion)){
      fimo$all_match[count2] = fimo$logs[count2]
      
      allMatches = allMatches + 1
      break;
    }
  }
  
}

topMatches <- 0
for (count1 in 1:length(chip_top$V1)) {
  for (count2 in 1:length(fimo$logs)) {
    if (chip_top$V1[count1] >= (fimo$start[count2] - rangeExtantion)
        && chip_top$V1[count1] <= (fimo$stop[count2] + rangeExtantion)){
      fimo$top_match[count2] = fimo$logs[count2]
      
      topMatches = topMatches + 1
      break;
    }
  }
  
}

#plot all results
plot(fimo$logs, type = 'l', main="FIMO with ChIP-seqs", 
     ylab = "-log of p-value", xlab = "index in FIMO", lwd=2)
points(fimo$top_match, col='red', pch=3, lwd=2)
points(fimo$all_match, col='blue', lwd=2)
legend(x='topright',legend = c('FIMO','TOP ChIP peaks','ALL ChIP peaks'), lty=c(1,0,0),pch = c(NA,3,1), col=c('black','red','blue'), lwd=c(2,2,2))

# text output
cat(topMatches, " out of " , length(chip_top$V1), " TOP ChIP regions matched with FIMO\n", sep = '')
cat(allMatches, " out of " , length(chip_all$V1), " ALL ChIP regions matched with FIMO", sep ='')



