setwd("/Users/sarah.killcoyne/Data/sky-cgh/output/26112012")
## Fragments ##
fg = read.table("fragments.txt", header=T, sep="\t")
fg = fg[order(fg$Chr),]

# Not using sex chromosomes generally, 
fg=fg[which(fg$Chr!="X" & fg$Chr!="Y"),]

## Lets ignore subbands (11.1) and just group them by the major band designation (11)
fg$Start = sub("\\.[0-9]+", "", fg$Start)
fg$End = sub("\\.[0-9]+", "", fg$End)

freq=table(fg$Chr)

# filter these out
for (i in 1:22)
	{
	
	print(i)
	print(nrow(fg$Chr==i))
	# p arm
	nt=fg[fg$Chr==i,]
	nt=nt[nt$Start!=ter[i, "p"] & nt$End!="pter",] 
	nt=nt[fg$Start!="pter" & fg$End!=ter[i, "p"],] 
#	fg=fg[which(fg$Chr==i & fg$Start!=ter[i, "p"] & fg$End!="pter"),] 
#	fg=fg[which(fg$Chr==i & fg$Start!="pter" & fg$End!=ter[i, "p"]),] 
	fg$Chr
	print(nrow(fg$Chr==i))
	}
