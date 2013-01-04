bp = read.table("breakpoints.txt", sep="\t", comment="#", header=T)

cancers = table(bp$Cancer)


write.table(cancers, row.names=F, col.names=F, quote=F, sep="\t")
write.table(cancers, file="~/Data/sky-cgh/output/cancers.txt", row.names=F, col.names=F, quote=F, sep="\t")