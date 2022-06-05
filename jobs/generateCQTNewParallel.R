library(batchtools)
library(snow)
library(tidyverse)


splitSaveRDSCqtCoeficientsFiles <- function(numberOfDevisions =10000) {
  rdsFiles <-  "data/cqtNew" %>%
    list.files( recursive = T,full.names = T, include.dirs = F)

  alreadyCalculatedFiles <- rdsFiles %>%
       lapply( function( file) file %>%
       readRDS() %>%
       select(path) ) %>%
       bind_rows() %>%
       rename(  filePath = path )


  listOfSets <- "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    read.csv() %>%
    mutate(filePath = filePath %>% as.character()) %>%
    anti_join(alreadyCalculatedFiles) %>%
    select(filePath) %>%
    unlist() %>%
    as.character() %>%
    split( c(1:(numberOfDevisions )))

  tibble(index = 1:numberOfDevisions, wavePath = listOfSets) %>%
    apply(  1,function( element) {

      spliterSaverBatchFunc <- function( wavePathes ) {
        library(tidyverse)
        library(reticulate)

        source_python( "cqtCalculation2.py")

        # cqtCalc = simpleCQT()

        generateCQTNpyFromFileAsList <- function(path) {
          path  %>%
            generateCQTNpyFromFile() %>%
            as.list() %>%
            set_names(c(1:13455) %>% as.character()) %>%
             c( list(path= path ))
        }
        wavePathes %>%  generateCQTNpyFromFileAsList() %>% bind_rows()

      }

      reg = makeRegistry(file.dir = NA, seed = 1)
      reg$cluster.functions = makeClusterFunctionsSocket(10)

      batchMap(fun = spliterSaverBatchFunc, element %>% magrittr::use_series(wavePath) %>% unlist() )
      submitJobs(resources = list(walltime = 3600, memory = 1024))
      getStatus()
      waitForJobs()

      fileIndex = (element %>% magrittr::use_series(index) ) + rdsFiles %>% length()
      reduceResults(function(result1, result2)  result1 %>% rbind(result2 ) ) %>%
        saveRDS(  sprintf( "data/cqtNew/cqtCoeficients_%s.rds",fileIndex )  )

    } )

}

splitSaveRDSCqtCoeficientsFiles(50000)