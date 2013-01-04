setwd("/Users/sarah.killcoyne/Data/sky-cgh/output/26112012")
## Ploidy changes ##
pdy = read.table("ploidy.txt", header=T, sep="\t")
pdy = pdy[order(pdy$Change),]
pdycnts = table(pdy$Change)
pdycnts = sort(pdycnts)
dev.new()
plot(pdycnts, type='h',  xaxt="n")
axis(1, at=1:length(pdycnts), lab=rownames(pdycnts))

#write.table(pdycnts[pdycnts > mean(pdycnts)],file="/Users/sarah.killcoyne/Data/sky-cgh/analysis/high-freq-ploidy.txt", quote=F, col.names=F)



