
setDataDirectory<-function(date = NA)
  {
  if(is.na(date))
    {
    date = format(Sys.Date(), format="%d%m%Y")
    }
  setwd(paste("~/Data/sky-cgh/output/", date, sep=""))
  message(getwd())
  }

resetwd<-function()
  {
  setwd("~/workspace/CancerCytogenetics")
  message(getwd())
  }