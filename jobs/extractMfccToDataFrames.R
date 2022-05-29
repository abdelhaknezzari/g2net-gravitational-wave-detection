library(batchtools)
library(snow)
library(tidyverse)


splitSaveCSVMfccFiles <- function(numberOfDevisions =10) {

  listOfSets <- "../machineLearningData/gravitationalWaves/mfcc/file_labels.csv"  %>%
    read.csv2() %>%
    nrow() %>%
    sample.int() %>%
    split( c(1:numberOfDevisions))



  tibble(index = 1:numberOfDevisions, listOfSets) %>%
    apply(  1,function( element) {

      spliterSaverBatchFunc <- function( indexes ) {
        source("functions.R")
        "../machineLearningData/gravitationalWaves/mfcc/file_labels.csv" %>%
          convertMFCCNpyFilesToDataFrameWith3Signals2(indexes)

      }

      reg = makeRegistry(file.dir = NA, seed = 1)
      reg$cluster.functions = makeClusterFunctionsSocket(10)

      batchMap(fun = spliterSaverBatchFunc, (element[['listOfSets']]  %>% split(  c(1:100))  ))
      submitJobs(resources = list(walltime = 3600, memory = 1024))
      getStatus()
      waitForJobs()

      reduceResults(function(result1, result2)  result1 %>% rbind(result2 ) ) %>%
      write.csv(  sprintf( "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.csv",element[['index']]) )



    } )



}

splitSaveCSVMfccFiles(200)

# testBatch <- function( startIndex ) {
#   source("functions.R")
#   "../machineLearningData/gravitationalWaves/mfcc/file_labels.csv" %>%
#     convertMFCCNpyFilesToDataFrameWith3Signals(startIndex  = startIndex, width =  10000)
# }
#
#
#
#
# "../machineLearningData/gravitationalWaves/mfcc/file_labels.csv" %>%
# convertMFCCNpyFilesToDataFrameWith3Signals2( indexes = listOfSets[[1]] ) %>%
# write.csv( "../machineLearningData/gravitationalWaves/mfcc/allMfccCoefficients1.csv")
#
#
# "../machineLearningData/gravitationalWaves/mfcc/file_labels.csv" %>%
#   convertMFCCNpyFilesToDataFrameWith3Signals2( indexes = listOfSets[[2]] ) %>%
#   write.csv( "../machineLearningData/gravitationalWaves/mfcc/allMfccCoefficients2.csv")
#






