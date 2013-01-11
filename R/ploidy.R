resetwd()
source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory()


## Ploidy changes ##
pdy = read.table("ploidy.txt", header=T, sep="\t")

sample_leukemia_bps = FALSE  # this makes no real difference in the ploidy frequencies
if (sample_leukemia_bps)
  {
  leuks = c('Acute myeloid leukemia', 'Acute lymphoblastic leukemia', "Non-hodgkin's lymphoma", 'Chronic myelogenous leukemia', 'Chronic lymphocytic leukemia')
  pdy = sampleCancers(pdy, "Cancer", leuks)
  }



colpdy = pdy[order(pdy$Ploidy),]
pdycnts = sort(table(pdy$Ploidy), decreasing=T)

plot(pdycnts, type='h',  xaxt="n")
axis(1, at=1:length(pdycnts), lab=rownames(pdycnts))

#write.table(pdycnts[pdycnts > mean(pdycnts)],file="/Users/sarah.killcoyne/Data/sky-cgh/analysis/high-freq-ploidy.txt", quote=F, col.names=F)



