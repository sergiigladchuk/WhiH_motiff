# all peak analysis
fimo <- read.csv('2016-12-05/top_pal/fimo.txt', sep='\t')
#View(fimo)

# palindrom correction
fimo <- fimo[-(seq(2,to=nrow(fimo),by=2)),]

# add ranges from all chip
chip_all <- read.csv('2016-12-05/all_ranges.txt', header = FALSE)

# add ranges from top chip
chip_top <- read.csv('2016-12-05/top_ranges.txt', header = FALSE)

#start vector
check <- c(rep('f',length(fimo$start)),rep('c',length(chip_top$V1)))
check <- as.data.frame(check)
check$pos <- c(fimo$start,chip_top$V1)

