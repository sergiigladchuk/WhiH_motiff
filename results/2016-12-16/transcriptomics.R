trans <- read.csv('2016-12-16/top10_pal_motif_counts.txt', sep='\t')
#str(trans)
#View(trans)

#get overall differance
trans$closest_match <- factor(trans$closest_match)

#split of ove and under expressed
#(-) indicates a decrease in expression of the gene in a whiH mutant compared to the wild-type; 
#(+) indicates an increase in expression of the gene in a whiH mutant compared to the wild-type.

#adding significant differences in transcriptomics +5% and -5% for each hour
signiff <- 0.05
times <- c ('8h', '10h', '12h', '14h', '16h', '18h', '20h')
outPut <- matrix(NA,nrow = 2, ncol = length(times))
outPutGrath <- matrix(NA,nrow = 8, ncol = length(times))
colnames(outPut) <- times
colnames(outPutGrath) <- times
rownames(outPut) <- c('OverExp. p-val.', 'UnderExp. p-val.')
rownames(outPutGrath) <- c('Sig. OverExp. Non-Motif', 'Sig. OverExp. Motif','Non-Sig. OverExp. Non-Motif', 'Non-Sig. OverExp. Motif', 
                           'Sig. UnderExp. Non-Motif','Sig. UnderExp. Motif', 'Non-Sig. UnderExp. Non-Motif','Non-Sig. UnderExp. Motif')
outPut <- as.data.frame(outPut)

#loop to goo through each time to get signifinace of motiff
col = 0
for (curentH in times) {
  col <- col + 1
  varName <- paste('affyLog_', curentH, sep='')
  trans.cur.over <- subset(trans, trans[[varName]] > 0)
  trans.cur.under <- subset(trans, trans[[varName]] < 0)
  trans.cur.over$sig <- ifelse(trans.cur.over[[varName]] >= quantile(trans.cur.over[[varName]],1 - signiff), 'S', 'N')
  trans.cur.under$sig <- ifelse(trans.cur.under[[varName]] <= quantile(trans.cur.under[[varName]],signiff), 'S', 'N')
  trans.cur.over$sig <- factor(trans.cur.over$sig)
  trans.cur.under$sig <- factor(trans.cur.under$sig)
  
  #View(trans.cur.over)
  #str(trans.cur.over)
  counts.over <- table(trans.cur.over$closest_match,trans.cur.over$sig)
  counts.under <- table(trans.cur.under$closest_match,trans.cur.under$sig)
  #plot(trans.8h.over$affyLog_8h~trans.8h.over$closest_match)
  
  #tests FIsher because can be lower than 4
  
  test.over <- fisher.test(counts.over, alternative = 'greater')
  test.under <- fisher.test(counts.under, alternative = 'greater')
  outPut[[curentH]][1] <- test.over$p.value
  outPut[[curentH]][2] <- test.under$p.value
  
  #for graph
  outPutGrath[1, col] <- counts.over[1,2]
  outPutGrath[2, col] <- counts.over[2,2]
  outPutGrath[3, col] <- counts.over[1,1]
  outPutGrath[4, col] <- counts.over[2,1]
  outPutGrath[5, col] <- counts.under[1,2]
  outPutGrath[6, col] <- counts.under[2,2]
  outPutGrath[7, col] <- counts.under[1,1]
  outPutGrath[8, col] <- counts.under[2,1]
  
}
cat('Statistical significance:')
outPut

cat('Counts for genes:')
outPutGrath
#graph
#outPutGrath <- t(outPutGrath)

#plot(outPutGrath[,1], type='b', xaxt='n', ylim=c(0,200), ylab=list('Counts',cex=1.5), xlab=list('Time',cex=1.5)
#     , col='blue', main=list('Count of significan over/under expresions',cex=1.5))
#axis(1, at=1:7, times)
#points(outPutGrath[,2], type='b', col='blue', pch = 5)
#points(outPutGrath[,3], type='b', col='red')
#points(outPutGrath[,4], type='b', col='red', pch = 5)

