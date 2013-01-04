# This script simply looks at the unparsed aberrations (e.g. der(12)t(12;14)(q15;q23)) and does some basic 
# frequency analysis.  Due to the inconsistent nature of aberration definition in karyotypes this is a very
# rough look at aberration frequency. What this show is a high frequency of ploidy aberrations (addition or deletion)
# of entire chromosomes as well as a few well characterized structural variations.

abr = read.table("aberrations.txt", header=T, sep="\t")
nrow(abr)
abr = abr[abr$Aberration != "r",]
abr = abr[abr$Aberration != "mar",]
nrow(abr)
abr_freq = table(abr$Aberration)
abr_freq = abr_freq[abr_freq>0]

range(abr_freq)
length(abr_freq)
sum(abr_freq[abr_freq==1])
## drop all of the singly occuring ones, they are interesting but at the moment just want to see the high frequency aberrations
abr_freq = abr_freq[abr_freq > 1]

std1 = sd(abr_freq)
high_freq = abr_freq[abr_freq>mean(abr_freq)+std1]

sorted_log = sort(log(high_freq), decreasing=T)

# take top 20 just to show
sorted_log[1:45]


plot(sorted_log[1:45], type="n")
text(sorted_log[1:45], labels=names(sorted_log[1:45]))

write.table(sort(high_freq, decreasing=T), quote=F, col.names=F, sep="\t", file="../../analysis/top-aberrations.txt")