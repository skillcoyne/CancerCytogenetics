
setDataDirectory<-function(date = NA)
  {
  if(is.na(date))
    {
    date = 'current'
    }
  setwd(paste("~/Data/sky-cgh/output/", date, sep=""))
  message(getwd())
  }

resetwd<-function()
  {
  setwd("~/workspace/CancerCytogenetics")
  message(getwd())
  }