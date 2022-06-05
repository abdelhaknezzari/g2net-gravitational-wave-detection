 keras_model_sequential() %>%
  layer_lstm(units = 63, return_sequences = TRUE, input_shape = c(19,32)) %>%
  layer_lstm(units = 32, return_sequences = TRUE) %>%
  layer_lstm(units = 32) %>% # return a single vector dimension 32
  layer_dense(units = 100) %>%
  layer_dropout(0.5) %>%
  layer_dense(1 ) %>%
  layer_activation("softmax") %>%
  compile(
    loss = 'categorical_crossentropy',
    optimizer = 'rmsprop',
    metrics = c('accuracy')
  ) %>% print()



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
  layer_dense(units = 1) %>%
  compile(loss = 'mae', optimizer = 'adam') %>%
  print()




 keras_model_sequential() %>%

   # Begin with 2D convolutional LSTM layer
   layer_conv_lstm_2d(
     input_shape = list(NULL,19,63,1),
     filters = 32, kernel_size = c(3,3),
     padding = "same",
     data_format="channels_last",
     return_sequences = TRUE
   ) %>%
   # Normalize the activations of the previous layer
   layer_batch_normalization() %>%

   # Add 3x hidden 2D convolutions LSTM layers, with
   # batch normalization layers between
   layer_conv_lstm_2d(
     filters = 40, kernel_size = c(3,3),
     padding = "same", return_sequences = TRUE
   ) %>%
   layer_batch_normalization() %>%
   layer_conv_lstm_2d(
     filters = 40, kernel_size = c(3,3),
     padding = "same", return_sequences = TRUE
   ) %>%
   layer_batch_normalization() %>%
   layer_conv_lstm_2d(
     filters = 40, kernel_size = c(3,3),
     padding = "same", return_sequences = TRUE
   ) %>%
   layer_batch_normalization() %>%

   # Add final 3D convolutional output layer
   layer_conv_3d(
     filters = 1, kernel_size = c(3,3,3),
     activation = "sigmoid",
     padding = "same", data_format ="channels_last"
   ) %>%
compile(
   loss = "binary_crossentropy",
   optimizer = "adadelta"
 ) %>% print()



 whitenWave <- function( wavePath, index) {
   sig1 =  wavePath %>% npyLoad() %>% .[index,]

   length1 = sig1 %>% length()

   spec1 = fft(sig1 *  hanning.w(length1 ) )

   mag1 = spec1 %>% Mod()

   (spec1 /mag1 ) %>%  ifft(  ) %>% Re() * sqrt( (length1 /2) )


 }



 library(Rwave)
library(signal)
 # https://pipiras.sites.oasis.unc.edu/timeseries/Nonstationary_2_-_Time_Frequency_-_Menu.html#what_is_this_all_about
 "../machineLearningData/gravitationalWaves/train/f/2/f/f2f0bbf138.npy" %>%
   whitenWave(2) %>%
#   npyLoad() %>% .[1,] %>%
   cwtp( noctave=16, nvoice=16) %>%
   image()
  plot()

sig1 =  "../machineLearningData/gravitationalWaves/train/f/2/f/f2f0bbf138.npy" %>%
   npyLoad() %>% .[1,]


 sig1 %>% specgram( n = min(1, 4096 ), Fs = 4096, window = hanning(4096),
          overlap = ceiling(length(window)/2))


 dir.exists
 dir.create


 "../machineLearningData/gravitationalWaves/train/ing_labels.csv" %>% read.csv() %>% View()


 source_python( "cqtCalculation2.py")


 "../machineLearningData/gravitationalWaves/train/f/f/f/fff0ac62d3.npy" %>%
   generateCQTNpyFromFile( ) %>%
   array_reshape(c(3,69,65)) %>% .[1,,] %>%
   image()

 "../machineLearningData/gravitationalWaves/train/f/f/f/fff0ac62d3.npy" %>%
   generateCQTNpyFromFile( ) %>%
   array_reshape(c(3,69,65)) %>% .[2,,] %>%
   image()

 "../machineLearningData/gravitationalWaves/train/f/f/f/fff0ac62d3.npy" %>%
   generateCQTNpyFromFile( ) %>%
   array_reshape(c(3,69,65)) %>% .[3,,] %>%  array_reshape(4485) %>% dim()
   image()


 generateCQTNpyFromFileAsList <- function(path =  "../machineLearningData/gravitationalWaves/train/f/f/f/fff0ac62d3.npy") {
    list( coeficients = path %>%
     generateCQTNpyFromFile() %>% toString(), wavePath = path ) }


 c("../machineLearningData/gravitationalWaves/train/f/f/f/fff1e83d2e.npy",
   "../machineLearningData/gravitationalWaves/train/f/f/f/fff0ac62d3.npy") %>%
   lapply( generateCQTNpyFromFileAsList ) %>%
   bind_rows() %>% View()




 "../machineLearningData/gravitationalWaves/train/f/f/f/fff0ac62d3.npy" %>%
   generateCQTNpyFromFile( ) %>% dim()
   array_reshape(c(3,69,65))


 "../machineLearningData/gravitationalWaves/train/f/f/f/fff1e83d2e.npy" %>%
   generateCQTNpyFromFile( ) %>%
   array_reshape(c(3,69,65)) %>% .[1,,] %>%
   image()

 "../machineLearningData/gravitationalWaves/train/f/f/f/fff1e83d2e.npy" %>%
   generateCQTNpyFromFile( ) %>%
   array_reshape(c(3,69,65)) %>% .[2,,] %>%
   image()

 "../machineLearningData/gravitationalWaves/train/f/f/f/fff1e83d2e.npy" %>%
   generateCQTNpyFromFile( ) %>%
   array_reshape(c(3,69,65)) %>% .[3,,] %>%
   image()




