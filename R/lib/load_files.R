
load.breakpoints<-function(files)
  {
  bp = NULL
  for (file in files)
    bp = rbind(bp, read.table(file, sep="\t", header=T))
  bp = bp[bp$total.karyotypes > 0,]

  bp$name = paste(bp$chr,bp$band,sep="")
  bp$name = sub("\\.[0-9]+","",bp$name)
  
  bp = bp[,c('chr','band','name','total.karyotypes')]
  
  for (x in names(which(table(bp$name) > 1)))
    {
    rows = which(bp$name == x)
    sum = sum(bp[rows, 'total.karyotypes'])
    newrow = bp[rows[grep("\\.[0-9]+",bp[rows, 'band'],invert=T)],][1,]
    newrow[,'total.karyotypes'] = sum
    bp = bp[-rows,]
    bp = rbind(bp,newrow)
    }
  bp = bp[order(-bp$total.karyotypes),] 
  return(bp)
  }


loadBands<-function(file = "../../genomic_info/chromosome_bands.csv")
  {
  t = read.table(file, sep=",", header=T)
  b = cbind(t$chromosome, t$band, t$start, t$end)
  b = b[order(b$chromosome),]  
  return(b)
  }
  
loadBreakpoints<-function(bpfile, use_sex = FALSE, ignore_subbands = TRUE)
  {
  bp = read.table(bpfile, sep="\t", comment="#", header=T)
  if (nrow(bp) <= 0) { stop("Breakpoint file is empty") }
  if (ignore_subbands)
    {
    ## Ignoring subbands (11.1) and just group them by the major band designation (11)
    #bp$Breakpoint = sub("\\.[0-9]+", "", bp$Breakpoint)
    bp$Breakpoint = clearSubbands(bp$Breakpoint)
    }
  if (!use_sex)
    {
    # Not using sex chromosomes generally, 
    bp = dropSexChr(bp)
    }
  return(bp)
  }



## if I care enough I can make one function for this stuff...
loadFragments<-function(file, use_sex = FALSE, ignore_subbands = TRUE)
  {
  ## Fragments ##
  fg = read.table(file, header=T, sep="\t")
  if (nrow(fg) <= 0) { stop("Fragments file is empty") }
  fg = fg[order(fg$Chr),]
  if (!use_sex)
    {
    # Not using sex chromosomes generally, 
    fg = dropSexChr(fg)
    }
  if (ignore_subbands)
    {
    ## Lets ignore subbands (11.1) and just group them by the major band designation (11)
    fg$Start = clearSubbands(fg$Start)
    fg$End = clearSubbands(fg$End)
    }
  return(fg)
  }

loadChromosomeInfo<-function(chrfile = "../../genomic_info/chromosome_gene_info_2012.txt")
  {
  chrinfo = read.table(chrfile, sep="\t", row.names=1, header=T)
  if (nrow(chrinfo) <= 0) { stop("Failed to load chromosome info file") }
  # don't need the mtDNA row
  chrinfo = chrinfo[ -(nrow(chrinfo)), ]
  
  for (i in 1:nrow(chrinfo))
    {
    row = chrinfo[i,]
    # leave pseudogenes out of it
    chrinfo[i,'Total.Prot.RNA'] = sum(row$Confirmed.proteins, row$Putative.proteins, row$miRNA, row$rRNA, row$snRNA, row$snoRNA, row$Misc.ncRNA)
    }
  return(chrinfo)
  }

dropSexChr<-function(df, colname = "Chr")
  {
  df = df[ which(df[[colname]]!="X" & df[[colname]]!="Y"),]
  return(df)
  }

clearSubbands<-function(col)
  {
  ## Ignoring subbands (11.1) and just group them by the major band designation (11)
  col = sub("\\.[0-9]+", "", col)
  return(col)
  }

sampleCancers<-function(df, colname = "Cancer", cancers)
  {
  message(nrow(df))
  dfcnc = df[ which(df[[colname]] %in% cancers ) ,]
  sampled = dfcnc[sample(1:nrow(dfcnc), 500, replace=FALSE),]
  noncnc = df[which(df[[colname]] %nin% cancers), ]
  df = rbind(noncnc, sampled)
  message(nrow(df))
  return(df)
  }
  


