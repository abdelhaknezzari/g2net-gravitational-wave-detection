dataGeneratorFromFiles3Waves <- function(files, labels,batchSize = 10 ) {
  fileIndex <- 0

  function() {
    if (fileIndex > length(files))
      fileIndex <<- 0

    nextIndex    <- fileIndex + batchSize
    rangeIndexes <- (fileIndex+1): nextIndex

    labelsToprocess <- labels[rangeIndexes ] %>% as.numeric()

    dataToProcess   <- files[ rangeIndexes ] %>%
      lapply(
        function(path) {
          wave1 <- path %>%
            readWaveFromNpy(1)  %>%
            calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
            keras::normalize(axis = -1, order = 2)
          wave2 <- path %>%
            readWaveFromNpy(2) %>%
            calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
            keras::normalize(axis = -1, order = 2)
          wave3 <- path %>%
            readWaveFromNpy(3) %>%
            calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
            keras::normalize(axis = -1, order = 2)
          list( wave1,wave2,wave3) %>% array_reshape( c(19,63,3) )
        } ) %>%
      array_reshape(c(batchSize,19,63,3))
    fileIndex <<- nextIndex
    list(dataToProcess, labelsToprocess )
  }
}



dataGeneratorFromFiles3WavesFromCSVs <- function( model ) {
  fileIndex <- 0

  function() {
    if (fileIndex >= 200 )
      fileIndex <<- 0

    nextIndex <- fileIndex + 1

    dataFromCsv <- "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.csv" %>%
      sprintf( nextIndex ) %>%
      read.csv()

    labelsToprocess <-  dataFromCsv  %>% select( target) %>% unlist() %>% as.numeric()
    dataToProcess   <-   dataFromCsv %>% select(-filePath,-id,-target,-X) %>% as.matrix() %>% array_reshape(c( ( dataFromCsv  %>% nrow() ) , 19,63,3))

    fileIndex <<- nextIndex
    list(dataToProcess, labelsToprocess )
  }
}



dataGeneratorFromFiles3WavesFromRdsTrain <- function(  ) {
  fileIndex <- 0

  function() {
    if (fileIndex >= 200 )
      fileIndex <<- 0

    nextIndex <- fileIndex + 1
    print(nextIndex)

    dataFromCsv <- "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.rds" %>%
      sprintf( nextIndex) %>%
      readRDS()

    numberOfEntries <- dataFromCsv %>% nrow()

    upperIndex = (numberOfEntries * 0.8 )  %>% as.integer()

    labelsToprocess <-  dataFromCsv  %>% select( target) %>% unlist() %>% as.numeric() %>% .[1:upperIndex]

    dataToProcess   <-   dataFromCsv %>% select(-filePath,-id,-target,-X) %>%
      slice( 1:upperIndex ) %>%
      as.matrix() %>% array_reshape(c( upperIndex , 19, 63, 3) )

    fileIndex <<- nextIndex
    list(dataToProcess, labelsToprocess )
  }
}


dataGeneratorFromFiles3WavesFromRdsValidate <- function(  ) {
  fileIndex <- 0

  function() {
    if (fileIndex >= 200 )
      fileIndex <<- 0

    nextIndex <- fileIndex + 1
    print(nextIndex)

    dataFromCsv <- "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.rds" %>%
      sprintf( nextIndex) %>%
      readRDS()

    numberOfEntries <- dataFromCsv %>% nrow()

    upperIndex = (numberOfEntries * 0.8 )  %>% as.integer()

    labelsToprocess <-  dataFromCsv  %>% select( target) %>%
      unlist() %>%
      as.numeric() %>%
      .[( upperIndex + 1):numberOfEntries ]


    dataSliced <-  dataFromCsv %>%
      select(-filePath,-id,-target,-X) %>%
      slice(( upperIndex + 1):numberOfEntries)

    dataToProcess   <-   dataSliced %>%
      as.matrix() %>%
      array_reshape(c(  dataSliced %>% nrow(), 19, 63, 3 ) )

    fileIndex <<- nextIndex
    list(dataToProcess, labelsToprocess )
  }
}






dataGeneratorFromFiles1Wave <- function(files, labels,batchSize = 10, waveIndex = 1 ) {
  fileIndex <- 0

  function() {
    if (fileIndex > length(files)) {
      fileIndex <<- 0
    }

    print("Hello")


    nextIndex = fileIndex + batchSize
    rangeIndexes <- (fileIndex+1): nextIndex
    labelsToprocess <- labels[ rangeIndexes ] %>% as.numeric()

    print(nextIndex)
    dataToProcess   <- files[ rangeIndexes ] %>%
      lapply(
        function(path) {

          wave <- path %>%
            readWaveFromNpy(waveIndex)  %>%
            calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
            keras::normalize(axis = -1, order = 2)
          wave %>% array_reshape( c(19,63) )
        } ) %>%
      array_reshape(c(batchSize,19,63))
    fileIndex <<- nextIndex
    list(dataToProcess, labelsToprocess )
  }
}


dataValidationGeneratorFromFiles1Wave <-  function (files, labels,batchSize = 10, waveIndex = 1 ) {
  fileIndex <- 0

  function() {
    if (fileIndex > length(files)) {
      fileIndex <<- 0
    }

    print("Hello Validation ")


    nextIndex = fileIndex + batchSize
    rangeIndexes <- (fileIndex+1): nextIndex
    labelsToprocess <- labels[ rangeIndexes ] %>% as.numeric()

    print(nextIndex)
    dataToProcess   <- files[ rangeIndexes ] %>%
      lapply(
        function(path) {

          wave <- path %>%
            readWaveFromNpy(waveIndex)  %>%
            calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
            keras::normalize(axis = -1, order = 2)
          wave %>% array_reshape( c(19,63) )
        } ) %>%
      array_reshape(c(batchSize,19,63))
    fileIndex <<- nextIndex
    list(dataToProcess, labelsToprocess )
  }
}


buildKerasModel3Waves <- function() {
  # Initialize sequential model
  model <- keras_model_sequential()

  model %>%

    # Start with hidden 2D convolutional layer being fed 32x32 pixel images
    layer_conv_2d(
      filter = 32, kernel_size = c(3,3), padding = "same",
      input_shape = c(19, 63, 3)
    ) %>%
    layer_activation("relu") %>%

    # Second hidden layer
    layer_conv_2d(filter = 32, kernel_size = c(3,3)) %>%
    layer_activation("relu") %>%

    # Use max pooling
    layer_max_pooling_2d(pool_size = c(2,2)) %>%
    layer_dropout(0.25) %>%

    # 2 additional hidden 2D convolutional layers
    layer_conv_2d(filter = 32, kernel_size = c(3,3), padding = "same") %>%
    layer_activation("relu") %>%
    layer_conv_2d(filter = 32, kernel_size = c(3,3)) %>%
    layer_activation("relu") %>%

    # Use max pooling once more
    layer_max_pooling_2d(pool_size = c(2,2)) %>%
    layer_dropout(0.25) %>%

    # Flatten max filtered output into feature vector
    # and feed into dense layer
    layer_flatten() %>%
    layer_dense(512) %>%
    layer_activation("relu") %>%
    layer_dropout(0.5) %>%

    # Outputs from dense layer are projected onto 10 unit output layer
    layer_dense(1 ) %>%
    layer_activation("softmax")

  opt <- optimizer_rmsprop(lr = 0.0001, decay = 1e-6)

  model %>% compile(
    loss = "categorical_crossentropy",
    optimizer = opt,
    metrics = "accuracy"
  )

  model %>% return()

}


buildKerasModel1Wave <- function() {
  # Initialize sequential model
  model <- keras_model_sequential()

  model %>%

    # Start with hidden 2D convolutional layer being fed 32x32 pixel images
    layer_conv_1d(
      filter = 32, kernel_size = 3, padding = "same",
      input_shape = c(19, 63)
    ) %>%
    layer_activation("relu") %>%

    # Second hidden layer
    layer_conv_1d(filter = 32, kernel_size = 3 ) %>%
    layer_activation("relu") %>%

    # Use max pooling
    layer_max_pooling_1d(pool_size = 2 ) %>%
    layer_dropout(0.25) %>%

    # 2 additional hidden 2D convolutional layers
    layer_conv_1d(filter = 32, kernel_size = 3, padding = "same") %>%
    layer_activation("relu") %>%
    layer_conv_1d(filter = 32, kernel_size = 3 ) %>%
    layer_activation("relu") %>%

    # Use max pooling once more
    layer_max_pooling_1d(pool_size = 2  ) %>%
    layer_dropout(0.25) %>%

    # Flatten max filtered output into feature vector
    # and feed into dense layer
    layer_flatten() %>%
    layer_dense(512) %>%
    layer_activation("relu") %>%
    layer_dropout(0.5) %>%

    # Outputs from dense layer are projected onto 10 unit output layer
    layer_dense(1 ) %>%
    layer_activation("softmax")

  opt <- optimizer_rmsprop(lr = 0.0001, decay = 1e-6)

  model %>% compile( loss = "categorical_crossentropy",   optimizer = opt,   metrics = "accuracy"   )
  #  model %>% compile( optimizer = "rmsprop",   loss = "categorical_crossentropy",   metrics = c("accuracy")  )

  model %>% return()

}



trainWith3Waves <- function(){

  labels <- readLabels()
  files <- readFiles()

  gen <- keras:::as_generator.function(dataGeneratorFromFiles3Waves(files, labels, batchSize = 10))

  model <- buildKerasModel3Waves()
  model %>%
    fit_generator(gen,
                  steps_per_epoch = ( files %>% length()) / 10, epochs = 200)
  model
}

trainWith3WavesWithMfccCSVs <- function(){
  model <- buildKerasModel3Waves()
  gen <- keras:::as_generator.function(dataGeneratorFromFiles3WavesFromCSVs( model))
  model %>%
    fit_generator(gen,
                  steps_per_epoch = 200 , epochs = 1 ,
                  verbose = 1)
  model
}

trainWith3WavesWithMfccRds <- function(){
  LossHistory <- R6::R6Class("LossHistory",
                             inherit = KerasCallback,
                             public = list(
                               losses = NULL,
                               on_batch_end = function(batch, logs = list()) {
                                 self$losses <- c(self$losses, logs[["loss"]])
                               }
                             ))
  # https://cran.r-project.org/web/packages/keras/vignettes/training_callbacks.html
  history <- LossHistory$new()

  model <- buildKerasModel3Waves()
  gen <- keras:::as_generator.function(dataGeneratorFromFiles3WavesFromRdsTrain( ))
  val <- keras:::as_generator.function(dataGeneratorFromFiles3WavesFromRdsValidate( ))

  model %>%
    fit_generator(gen,
                  steps_per_epoch = 200 , epochs = 1,
                  callbacks = list(callback_csv_logger("log.csv", separator = ",", append = FALSE)),
                  # callback_csv_logger("log.csv", separator = ",", append = FALSE)
                  validation_data = val ,
                  validation_steps = 200
    )
}


splitDataTrainValidate <- function() {
  files <- "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves()


  numberRows <- files %>% nrow()
  indexes    <- sample( numberRows , replace = T )
  upperTrainingIndex <- (numberRows*0.8)

  trainIndexes <- indexes[ 1: upperTrainingIndex ]
  validIndexes <- indexes[ upperTrainingIndex: numberRows ]

  trainData <- files[trainIndexes,] %>%
    select(filePath) %>%
    unlist() %>%
    as.character()

  trainLabels <- files[trainIndexes,] %>%
    select(target) %>%
    unlist() %>%
    as.character()



  validData <- files[validIndexes,] %>%
    select(filePath) %>%
    unlist() %>%
    as.character()

  validLabels <- files[validIndexes,] %>%
    select(target) %>%
    unlist() %>%
    as.character()

  list(trainData, trainLabels , validData , validLabels  )



}

trainWith1Wave <- function(waveIndex = 1, batchSize =200){

  splitData <- splitDataTrainValidate()


  gen <- keras:::as_generator.function(dataGeneratorFromFiles1Wave(splitData[[1]], splitData[[2]], batchSize = batchSize, waveIndex = waveIndex))
  val <- keras:::as_generator.function(dataValidationGeneratorFromFiles1Wave(splitData[[3]],splitData[[4]], batchSize = batchSize, waveIndex = waveIndex))

  stepsTrain <- ( splitData[[1]] %>% length() / batchSize ) %>% as.integer()
  stepsValid <- ( splitData[[3]] %>% length() / batchSize ) %>% as.integer()

  model <- buildKerasModel1Wave()
  model %>%
    fit_generator( gen,
                   epochs = 1,
                   steps_per_epoch = stepsTrain,
                   callbacks = list(callback_progbar_logger(count_mode = "samples"),
                                    #callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.1) ,
                                    callback_remote_monitor(                           root = "https://localhost:9000",
                                                                                       path = "/publish/epoch/end/",
                                                                                       field = "data",
                                                                                       headers = NULL,
                                                                                       send_as_json = FALSE )
                   ),
                   validation_data = val ,
                   validation_steps = stepsValid)
  model

}



model2With3Waves <- function() {
  keras_model_sequential() %>%
    layer_lstm(units = 3,
               input_shape = c(63, 3000),
               batch_size = 19,
               return_sequences = TRUE,
               stateful = TRUE) %>%
    layer_dropout(rate = 0.5) %>%
    layer_lstm(units = 50,
               return_sequences = FALSE,
               stateful = TRUE) %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(1 ) %>%
    layer_activation("softmax") %>%
    compile(loss = 'mae', optimizer = 'adam')

}

trainWith3WavesWithMfccRds <- function() {
  LossHistory <- R6::R6Class("LossHistory",
                             inherit = KerasCallback,
                             public = list(
                               losses = NULL,
                               on_batch_end = function(batch, logs = list()) {
                                 self$losses <- c(self$losses, logs[["loss"]])
                               }
                             ))
  # https://cran.r-project.org/web/packages/keras/vignettes/training_callbacks.html
  history <- LossHistory$new()

  model <- model2With3Waves()
  gen <- keras:::as_generator.function(dataGeneratorFromFiles3WavesFromRdsTrain( ))
  val <- keras:::as_generator.function(dataGeneratorFromFiles3WavesFromRdsValidate( ))

  model %>%
    fit_generator(gen,
                  steps_per_epoch = 200 , epochs = 1,
                  callbacks = list(callback_csv_logger("log.csv", separator = ",", append = FALSE)),
                  # callback_csv_logger("log.csv", separator = ",", append = FALSE)
                  validation_data = val ,
                  validation_steps = 200
    )
}


trainValidateByOneWave <- function( dimIndex) {


  dataGeneratorFromFiles1WavesFromRdsTrain <- function( dimIndex ) {
    fileIndex <- 0

    function() {
      if (fileIndex >= 200 )
        fileIndex <<- 0

      nextIndex <- fileIndex + 1
      print(nextIndex)

      dataFromCsv <- "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.rds" %>%
        sprintf( nextIndex) %>%
        readRDS()

      numberOfEntries <- dataFromCsv %>% nrow()

      upperIndex = (numberOfEntries * 0.8 )  %>% as.integer()

      labelsToprocess <-  dataFromCsv  %>% select( target) %>% unlist() %>% as.numeric() %>% .[1:upperIndex]

      dataToProcess   <-   dataFromCsv %>% select(-filePath,-id,-target,-X) %>%
        slice( 1:upperIndex ) %>%
        as.matrix() %>% array_reshape(c( upperIndex , 19, 63, 3) ) %>%
       .[,,,dimIndex]

      fileIndex <<- nextIndex
      list(dataToProcess, labelsToprocess )
    }
  }


  dataGeneratorFromFiles1WavesFromRdsValidate <- function( dimIndex ) {
    fileIndex <- 0

    function() {
      if (fileIndex >= 200 )
        fileIndex <<- 0

      nextIndex <- fileIndex + 1
      print(nextIndex)

      dataFromCsv <- "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.rds" %>%
        sprintf( nextIndex) %>%
        readRDS()

      numberOfEntries <- dataFromCsv %>% nrow()

      upperIndex = (numberOfEntries * 0.8 )  %>% as.integer()

      labelsToprocess <-  dataFromCsv  %>% select( target) %>%
        unlist() %>%
        as.numeric() %>%
        .[( upperIndex + 1):numberOfEntries ]


      dataSliced <-  dataFromCsv %>%
        select(-filePath,-id,-target,-X) %>%
        slice(( upperIndex + 1):numberOfEntries)

      dataToProcess   <-   dataSliced %>%
        as.matrix() %>%
        array_reshape(c(  dataSliced %>% nrow(), 19, 63, 3 ) ) %>%
        .[,,,dimIndex]

      fileIndex <<- nextIndex
      list(dataToProcess, labelsToprocess )
    }
  }

  getLstmModel <- function() {
    keras_model_sequential() %>%
      layer_lstm(units = 32, return_sequences = TRUE, input_shape = c(19,63)) %>%
      layer_lstm(units = 63, return_sequences = TRUE) %>%
      layer_lstm(units = 63) %>% # return a single vector dimension 32
      layer_dense(units = 100) %>%
      layer_dropout(0.5) %>%
      layer_dense(1 ) %>%
      layer_activation("softmax") %>%
      compile(
        loss = 'categorical_crossentropy',
        optimizer = 'rmsprop',
        metrics = c('accuracy')
      )
  }

  modelForOneWave <- getLstmModel()

  gen <- keras:::as_generator.function(dataGeneratorFromFiles1WavesFromRdsTrain(dimIndex))
  val <- keras:::as_generator.function(dataGeneratorFromFiles1WavesFromRdsValidate(dimIndex))

  modelForOneWave  %>%
    fit_generator( gen,
                   epochs = 4,
                   callbacks = list(callback_csv_logger("log.csv", separator = ",", append = FALSE)),
                   steps_per_epoch = 200,
                   validation_data = val ,
                   validation_steps = 200)


}

