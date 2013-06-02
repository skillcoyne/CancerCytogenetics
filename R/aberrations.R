# This script simply looks at the unparsed aberrations (e.g. der(12)t(12;14)(q15;q23)) and does some basic 
# frequency analysis.  Due to the inconsistent nature of aberration definition in karyotypes this is a very
# rough look at aberration frequency. What this show is a high frequency of ploidy aberrations (addition or deletion)
# of entire chromosomes as well as a few well characterized structural variations.
resetwd()
source("lib/wd.R")
source("lib/load_files.R")
`%nin%`=Negate(`%in%`) 
setDataDirectory()

datadir = "~/Data/sky-cgh/output/current"
setwd(datadir)

#total_samples = totalKaryotypes()

df = read.table("aberrations.txt", header=T, sep="\t")
nrow(df)

# Not dealing with these aberrations currently
`%nin%`=Negate(`%in%`)
df = df[df$class %nin% c('ring', 'mar', 'loss', 'gain'),]

leuk_rows = which(df$leukemias > 0)
non_leuk = which(df$leukemias <= 0)

sample_leukemia_bps = FALSE
samples = df

if (sample_leukemia_bps)
  {
  message("Sampling leukemia karyotypes")
  leukemias = sample(leuk_rows, 500, replace=F)
  samples = rbind( df[leukemias,], df[non_leuk,]  )
  }
nrow(df)
nrow(samples)

# grab all aberrations with more than 1 karyotype
abrs = samples[,c(2,4)]
abrs = abrs[abrs[,2] > 1,]

# turn it into a table
abr_freq = vector(mode="numeric", length=nrow(abrs))
names(abr_freq) = abrs[,1]
abr_freq[1:length(abr_freq)] = abrs[,2]

range(abr_freq)
length(abr_freq)

## drop all of the singly occuring ones, they are interesting but at the moment just want to see the high frequency aberrations
#abr_freq = abr_freq[abr_freq > 1]

std1 = sd(abr_freq)
high_freq = abr_freq[abr_freq>mean(abr_freq)+std1]

filename = "aberration_samples_all.txt"
if (sample_leukemia_bps) filename = "aberration_freq_sampled.txt"

write.table(sort(high_freq, decreasing = T), col.names=F, sep="\t", quote=F, file=filename)

sorted_log = sort(log(high_freq), decreasing=T)

# take top 20 just to show
sorted_log[1:45]


plot(sorted_log[1:45], type="h")
text(sorted_log[1:45], labels=names(sorted_log[1:45]))

