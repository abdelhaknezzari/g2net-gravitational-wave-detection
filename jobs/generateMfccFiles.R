source('functions.R')
readFiles() %>% lapply( function( filePath) {
  filePath %>% save3SignalMFCCAsNpy()
} )




