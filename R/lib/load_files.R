totalKaryotypes<-function()
  {
  t = read.table("totals.txt", sep="\t", row.names=1, header=T)
  return( sum(t$Total) )
  }

unknowns<-function()
  {
  t = read.table("totals.txt", sep="\t", row.names=1, header=T)
  return( sum(t$Unknown) )
  }

loadBreakpoints<-function(bpfile, use_sex = FALSE, ignore_subbands = TRUE)
  {
  bp = read.table(bpfile, sep="\t", comment="#", header=T)
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