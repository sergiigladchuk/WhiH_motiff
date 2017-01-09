pref <- 'best_top_discrim_motif'
trans <- read.csv(paste('2017-01-03/',pref,'_genes.TAB',sep = ''), sep='\t')
#generate link to genes table
#kable(paste('[Link to whole match table of genes](2017-01-03/',pref,'_genes.TAB)',sep = ''))
#get overall differance
trans$closest_match <- factor(trans$closest_match)

#to clear the file for first run
#write('Files',file = '2017-01-03/rank_over.txt',append = FALSE);
#write('Files',file = '2017-01-03/rank_under.txt',append = FALSE);
#split of ove and under expressed
#(-) indicates a decrease in expression of the gene in a whiH mutant compared to the wild-type; 
#(+) indicates an increase in expression of the gene in a whiH mutant compared to the wild-type.

#adding significant differences in transcriptomics +5% and -5% for each hour
times <- c ('8h', '10h', '12h', '14h', '16h', '18h', '20h')
outPut <- matrix(NA,nrow = 4, ncol = length(times))
colnames(outPut) <- times
rownames(outPut) <- c('OverExp. flagged genes','OverExp. p-val.', 'UnderExp. flagged genes','UnderExp. p-val.');

#loop to goo through each time to get signifinace of motiff
colCount <- 0

for (curentH in times) {
  colCount <- colCount + 1
  
  varName <- paste('affyLog_', curentH, sep='')
  trans.cur.over <- subset(trans, trans[[varName]] > 0)
  trans.cur.under <- subset(trans, trans[[varName]] < 0)
  
  #write tables with ouput for flaged genes and generate links
  fnameOver <- paste('2017-01-03/',pref,'_lists/over_',curentH,'.TAB',sep = '')
  fnameUnder <- paste('2017-01-03/',pref,'_lists/under_',curentH,'.TAB',sep = '')
  
  #outframes
  outOver <- cbind(trans.cur.over[c(1:10)],trans.cur.over[[varName]])
  outOver <- subset(outOver,outOver$closest_match == '1')
  outUnder <- cbind(trans.cur.under[c(1:10)],trans.cur.under[[varName]])
  outUnder <- subset(outUnder,outUnder$closest_match == '1')
  #rename columns
  colnames(outOver) <- c(colnames(trans)[1:10],varName)
  colnames(outUnder) <- c(colnames(trans)[1:10],varName)
  #sort
  outOver <- outOver[rev(order(outOver[[varName]])),]
  outUnder <- outUnder[order(outUnder[[varName]]),]
  #output
  write.table(outOver, file = fnameOver,sep = '\t', row.names = FALSE, na='')
  write.table(outUnder, file = fnameUnder,sep = '\t', row.names = FALSE, na='')
  
  
  #Mann‐Whitney non-parametric independent samples test
  mann.whitney.over <- wilcox.test(trans.cur.over[[varName]]~trans.cur.over$closest_match, alternative='less')
  mann.whitney.under <- wilcox.test(trans.cur.under[[varName]]~trans.cur.under$closest_match, alternative='greater')
  outPut[1,colCount] <- paste0('[',length(which(trans.cur.over$closest_match == '1')),' out of ',length(trans.cur.over$closest_match),'](',fnameOver,')')
  outPut[2,colCount] <- format(mann.whitney.over$p.value, scientific = -2, digits = 4)
  #add color to significant p-val
  if (mann.whitney.over$p.value < 0.05) {
    outPut[2,colCount] <- paste0('<span style="background-color:lightgreen">',outPut[2,colCount],'</span>')
    #add list name of significant genes to total ranking list
    write(paste0('results/',fnameOver),file = '2017-01-03/rank_over.txt',append = TRUE);
  }
  outPut[3,colCount] <- paste0('[',length(which(trans.cur.under$closest_match == '1')),' out of ',length(trans.cur.under$closest_match),'](',fnameUnder,')')
  outPut[4,colCount] <- format(mann.whitney.under$p.value, scientific = -2, digits = 4)
  if (mann.whitney.under$p.value < 0.05) {
    outPut[4,colCount] <- paste0('<span style="background-color:lightgreen">',outPut[4,colCount],'</span>')
    write(paste0('results/',fnameUnder),file = '2017-01-03/rank_under.txt',append = TRUE);
  }
 
}

# library(Gmisc)
# htmlTable(outPut)
# kable(outPut)
