source('functions.R')
library(reticulate)
source_python( "cqtCalculation2.py")

readFiles() %>% lapply( function( filePath) {
  filePath %>% save3SignalCQTAsNpy()
} )




