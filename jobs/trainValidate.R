source('functions.R')
library(batchtools)
library(tidyverse)


splitDataTrainValidate <- function( data, job, index, ...) {
  library(tidyverse)

  mfccData <- sprintf( "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.rds",index) %>%
    readRDS()
  sizeOfData <- mfccData %>%
    nrow()
  indexes <- sizeOfData %>%
    sample.int()
  upperIndex <- sizeOfData * 0.8 %>% as.numeric()
  list( train = mfccData %>% slice(  indexes[1:upperIndex] ), test = mfccData %>% slice(  indexes[upperIndex:sizeOfData] )  )
}

svm.wrapper = function(data, job, instance, ...) {
  library("e1071")
  mod = svm(target ~ ., data = instance$train,  kernel = "linear")
  pred = predict(mod, newdata = instance$test, type = "class")
  table(instance$test$target, pred)
}


forest.wrapper = function(data, job, instance, ...) {
  library(ranger)
  mod = ranger(target ~ ., data = instance$train, write.forest = TRUE)
  pred = predict(mod, data = instance$test )
  table(instance$test$target, pred$predictions)
}

reg = makeExperimentRegistry(file.dir = NA, seed = 1)
reg$cluster.functions = makeClusterFunctionsSocket(10)

batchtools::addProblem(name = "mfccData", fun = splitDataTrainValidate, seed = 42)
batchtools::addAlgorithm(name = "svm", fun = svm.wrapper)
batchtools::addAlgorithm(name = "forest", fun = forest.wrapper)

library(data.table)
pdes = list(mfccData = data.table(index = 1:10 ))

ades = list(
  svm = CJ(kernel = c( "linear" ), epsilon = c(0.05)),
  forest = data.table(ntree = c(500))
)

addExperiments(pdes, ades, repls = 5)

submitJobs()

waitForJobs()


model <- buildKerasModel3Waves()



vv <- 1:20 %>%  lapply( function(index)  {
  sprintf( "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.rds",index) %>% readRDS()
}  )  %>%
  bind_rows()


instance <- list( train = vv[1: 44800,] , test=  vv[ 44801:56000, ]  )



vv <- 1:20 %>% lapply( function( fileIndex) {

  trainvalData <- splitDataTrainValidate(fileIndex, 0.8)

  history <- model %>% fit(
    trainvalData$train %>% select(-filePath,-id,-target,-X) %>% as.matrix() %>% array_reshape(c( (  trainvalData$train  %>% nrow() ) , 19,63,3)),
    trainvalData$train %>% select( target) %>% unlist() %>% as.numeric(),
    # validation_data = list(
    #   trainvalData$validate %>% select(-filePath,-id,-target,-X) %>% as.matrix() %>% array_reshape(c(( trainvalData$validate %>% nrow() ) , 19,63,3)),
    #   trainvalData$validate %>% select( target) %>% unlist() %>% as.numeric()
    # ),
    callbacks = list(
      callback_early_stopping(patience = 5),
      callback_reduce_lr_on_plateau(patience = 3)
    ),

    batch_size = 200,
    epochs = 140
  )

  history
})




submitJobs(resources = list(walltime = 3600, memory = 1024))
getStatus()
waitForJobs()