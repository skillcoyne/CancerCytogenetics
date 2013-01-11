source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory(date = '09012013')


bp = loadBreakpoints("breakpoints.txt")

cnc = read.table("cancers.txt", sep="\t", header=T)


ncbi = cnc[ which(cnc$Source == 'ncbi'), ]
cnt = table(ncbi$Cancer)


sc = table(cnc$Source, cnc$Cancer)

barplot(sc, col=c("darkblue", "red", "green"))