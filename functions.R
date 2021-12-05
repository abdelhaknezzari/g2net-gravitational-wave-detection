library(RcppCNPy)
library(tidyverse)
library(seewave)
library(tuneR)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)
library(scatterplot3d)
library(Rtsne)
library(keras)



readWaveFromNpy <- function(filePath,index) {
  filePath %>%
    npyLoad() %>%
    .[index,] %>%
    scale()  %>%
    Wave( left = .,samp.rate = 2048 , bit = 32 )
}

# x label formatter
s_formatter <- function(x){
  lab <- paste0(x, " s")
}

# y label formatter
khz_formatter <- function(y){
  lab <- paste0(y, " kHz")
}

buildDataFrameOfWave <- function(wav) {
  sample <- seq(1:length(wav@left))
  time <- sample/wav@samp.rate
  sample.left <- wav@left %>% cbind() %>% as.vector()
  sample %>% data.frame( time, sample.left)
}


buildSpectrograph <- function(wav,fmax = 0.06,wl = 2048,wn  = "hamming",ovlp = 0) {
  hot_colors <- inferno(n=20)

  hot_theme_grid <- theme(panel.grid.major.y = element_line(color="black", linetype = "dotted"),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(fill="transparent"),
                          panel.border = element_rect(linetype = "solid", fill = NA, color = "grey"),

                          axis.line = element_blank(),

                          # legend.position = "top",
                          # legend.justification = "right",
                          # legend.background = element_rect(fill="black"),
                          # legend.key.width = unit(50, "native"),
                          # legend.title = element_text(size=16, color="grey"),
                          # legend.text = element_text(size=16, color="grey"),

                          plot.background = element_rect(fill="black"),
                          plot.margin = margin(1,1,0,1, "lines"),

                          axis.title = element_blank(),
                          axis.text = element_text(size=16, color = "grey"),
                          axis.text.x = element_blank(),
                          axis.ticks = element_line(color="grey"))

  wav %>%
    ggspectro( flim=c(0,fmax), wl = wl, wn = wn,zp=0 , ovlp = ovlp )+
    scale_x_continuous(labels=s_formatter, expand = c(0,0))+
    scale_y_continuous(breaks = seq(from = 5, to = 20, by=5), expand = c(0,0), labels = khz_formatter, position = "right")+
    geom_raster(aes(fill=amplitude), hjust = 0, vjust = 0, interpolate = F)+
    scale_fill_gradientn(colours = hot_colors, name = "Amplitude \n (dB)", na.value = "transparent", limits = c(-60,0))+
    hot_theme_grid

}

buildOscilograph <- function(df) {
  oscillo_theme_dark <- theme(panel.grid.major.y = element_line(color="black", linetype = "dotted"),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_rect(fill="transparent"),
                              panel.border = element_rect(linetype = "solid", fill = NA, color = "grey"),
                              axis.line = element_blank(),
                              legend.position = "none",
                              plot.background = element_rect(fill="black"),
                              plot.margin = unit(c(0,1,1,1), "lines"),
                              axis.title = element_blank(),
                              axis.text = element_text(size=14, color = "grey"),
                              axis.ticks = element_line(color="grey"))
  df %>%
    ggplot() +
    geom_line(mapping = aes(x=time, y=sample.left), color="grey")+
    scale_x_continuous(labels=s_formatter, expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), position = "right")+
    geom_hline(yintercept = 0, color="white", linetype = "dotted")+
    oscillo_theme_dark
}


buildFFT <- function( wav ,fmax = 0.06,wl = 2048,wn = "hamming") {

  oscillo_theme_dark <- theme(panel.grid.major.y = element_line(color="black", linetype = "dotted"),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_rect(fill="transparent"),
                              panel.border = element_rect(linetype = "solid", fill = NA, color = "grey"),
                              axis.line = element_blank(),
                              legend.position = "none",
                              plot.background = element_rect(fill="black"),
                              plot.margin = unit(c(0,0,0,0), "lines"),
                              axis.title = element_blank(),
                              axis.text = element_text(size=14, color = "grey"),
                              axis.ticks = element_line(color="grey"))


  wav %>%
    spec(  plot = 0, wl= wl,wn = wn ) %>%
    as.data.frame() %>%
    filter( x <= fmax ) %>%
    ggplot() +
    geom_line(mapping = aes(x=x, y=y), color="grey")+
    scale_x_continuous(labels=khz_formatter, expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), position = "right")+
    geom_hline(yintercept = 0, color="white", linetype = "dotted")+
    coord_flip()+
    oscillo_theme_dark
}




#-------------------------------------------
## PLOT GRID
#-------------------------------------------

showPlotsInGrid <- function( hotplot , oscplot , fftPlot ) {
  gA <- hotplot %>% ggplot_build() %>% ggplot_gtable()
  gB <- oscplot %>% ggplot_build() %>% ggplot_gtable()
  gC <- fftPlot %>% ggplot_build() %>% ggplot_gtable()

  maxWidth = grid::unit.pmax(gA$widths, gB$widths)

  gA$widths <- as.list(maxWidth)
  gB$widths <- as.list(maxWidth)

  layo <- rbind(c(1,1,3),
                c(1,1,3),
                c(1,1,3),
                c(2,2,NA))

  grid.newpage()
  grid.arrange(gA, gB,gC, layout_matrix = layo)
}


showPlots <- function(wav,fmax = 0.06,wl = 2048,wn = "hamming",ovlp =0 ) {
  hotplot <- wav %>% buildSpectrograph(fmax = fmax,wl = wl,wn = wn,ovlp = ovlp)
  oscplot <- wav %>% buildDataFrameOfWave() %>% buildOscilograph()
  fftPlot <- wav %>% buildFFT(fmax = fmax,wl = wl,wn = wn)
  showPlotsInGrid( hotplot, oscplot,fftPlot)
}



spectreTest <- function(fname) {
  fname %>% npyLoad() %>% scale()-> imat
  imat <- imat * 20
  imat %>%  .[2,] %>% Wave( left = .,samp.rate = 2048 )  %>%
    spectro( flim = c(0,0.2) , palette = colorRampPalette(c("white", "blue", "green")) , ovlp = 50)

}

waveAuditorySpectrumScale <- function( wave,wl = 512) {
  # sampling frequency
  fs <- wave@samp.rate
  #  Power spectrum using hamming window
  wave@left %>% powspec( sr=fs,wintime=wl/fs, steptime=wl/fs) %>%
    #     frequency band conversion -> reduce frequencies to auditory frequency scale
    audspec(sr=fs, minfreq=0, maxfreq=fs/2,   nfilts=26, fbtype="htkmel")

}

scaleAndPlotAudSpecToMelFrequencies <- function(wave) {
  audspecWave <-  wave %>% waveAuditorySpectrumScale()
  fs <- wave@samp.rate
  # time scale
  at <- seq(0, 1, length=5)
  time <- round(seq(0, duration(wave), length=5), 1)
  # Hz frequency scale
  hz <- round(seq(fs/512, fs/2, length=5))
  # mel frequency scale
  mel <- hz %>%  hz2mel( htk=TRUE) %>% round()
  # plot
  par(mar=c(5.1, 4.1, 4.1, 4.1), las=1)
  col <- gray((512:0)/512)
  audspecWave$aspectrum %>% t() %>% image( col=col,
                                           axes=FALSE, xlab="Time (s)", ylab="Frequency (mel)")
  axis(side=1, at=at, labels=time)
  axis(side=2, at=at, labels=mel)
  axis(side=4, at=0:25/25, labels=1:26,)
  mtext("Mel-frequency filter #", side=4, las=0, line=2.5)
  abline(h=(0:25/25)+1/(25*2), col="lightgray")
  abline(v=(0:73/73)+1/(73*2), col="lightgray")
  box()
}

cepstralCoefs <- function(wave) {
  audspecWave <-  wave %>%
    waveAuditorySpectrumScale() %>%
    .$aspectrum %>%
    spec2cep( ncep=13, type="t3")
}

scaleAndPlotCepstralCoefs <- function(wave) {
  cepstralCoef <-  wave %>% cepstralCoefs()
  fs <- wave@samp.rate
  # time scale
  at <- seq(0, 1, length=5)
  time <- round(seq(0, duration(wave), length=5), 1)
  # Hz frequency scale
  hz <- round(seq(fs/512, fs/2, length=5))
  # mel frequency scale
  mel <- hz %>%  hz2mel( htk=TRUE) %>% round()
  # plot
  par(mar=c(5.1, 4.1, 4.1, 4.1), las=1)
  col <- gray((512:0)/512)
  cepstralCoef$cep %>%
    t() %>%
    image( col=col,axes=FALSE, xlab="Time (s)", ylab="Frequency (mel)")
  axis(side=1, at=at, labels=time)
  axis(side=2, at=at, labels=mel)
  axis(side=4, at=0:25/25, labels=1:26,)
  mtext("Mel-frequency filter #", side=4, las=0, line=2.5)
  abline(h=(0:25/25)+1/(25*2), col="lightgray")
  abline(v=(0:73/73)+1/(73*2), col="lightgray")
  box()
}


mfccCoefs <- function(wave) {
  fs <- wave@samp.rate # sampling frequency
  wl <- 512 # STDFT window size
  ncep <- 13 # final number of MFCCs
  wave %>%
    preemphasis( alpha=0.97, output="Wave") %>%
    .@left %>%
    powspec( sr=fs, wintime=wl/fs, steptime=wl/fs) %>%
    audspec( sr=fs, nfilts=ncep*2, fbtype="htkmel") %>%
    .$aspectrum %>%
    spec2cep( ncep=ncep, type="t3") %>%
    .$cep %>%
    lifter( lift=ncep-1, htk=TRUE)
}

mfccCoefs2 <- function(wave,ncep = 13,wl = 512, fbtype="htkmel" , dcttype="t3") {
  fs <- wave@samp.rate # sampling frequency
  wave %>%  melfcc( sr=fs,
                    wintime=wl/fs, hoptime=wl/fs,
                    numcep=ncep, nbands=ncep*2,
                    fbtype= fbtype, dcttype=dcttype,
                    htklifter=TRUE, lifterexp=ncep-1,
                    frames_in_rows=FALSE,
                    spec_out=TRUE)
}

scalePlotMfcc <- function(wave) {
  ## time scale
  at <- seq(0, 1, length=5)
  time <- wave %>% duration() %>% seq(0, ., length=5) %>% round(1)
  ## plot
  # grey scale of colors
  col <- ((512:0)/512) %>% gray()
  par(las=1)
  wave %>%
    mfccCoefs2() %>%
    .$cepstra %>%
    t() %>%
    image(col=col,axes=FALSE, xlab="Time (s)", ylab="MFCC #")
  axis(side=1, at=at, labels=time)
  axis(side=2, at=0:12/12, labels=1:13)
  abline(h=(0:12/12)+1/(12*2), col="lightgray")
  abline(v=(0:73/73)+1/(73*2), col="lightgray")
  box()
}

getWaveTitle <- function( waveIndex ) {

  switch( waveIndex,
     "LIGO Hanford" ,
    "LIGO Livingston",
    "Virgo"
  )
}

ggplotMFCC <- function(wave, waveIndex = 1, waveClass = 1 , filePath = "", ncep = 13,wl = 512, fbtype="htkmel", dcttype="t3" ) {
  signalDuration <- wave %>% duration()
  filName <- filePath %>% substr(start = 58, stop = 71 )
  wave %>%
    mfccCoefs2(ncep,wl, fbtype, dcttype) %>%
    .$cepstra %>%
    reshape2::melt() %>%
    mutate( Var2 = signalDuration * (Var2/max(Var2))  ) %>%
    ggplot( aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    scale_fill_gradientn(colours = ((512:0)/512) %>% gray()  ) +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks=13) +
    labs(x="Time (s)", y="Coeficients", title= sprintf( "MFCC coeficients: %s, Class: %s, Path: %s", getWaveTitle(waveIndex), waveClass, filName  ) ) +
    theme_bw() +
    theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
          axis.text.y=element_text(size=9),
          plot.title=element_text(size=11))
}


getListOfWaves <- function(filePath )  {
  filePath %>%
    read.csv() %>%
    as_tibble()
}


getRandomWaveFilePath <- function ( filePath , dataClass , parentPath = "") {
  filePath %>%
    getListOfWaves() %>%
    filter(target == dataClass  ) %>%
    sample_n( size = 10 ) %>%
    select(-X) %>%
    select(filePath) %>%
    slice(1) %>%
    unlist() %>%
    as.character() %>%
    paste0(parentPath,.)
}


plotOneWaveMFCC <- function(filePath, dataClass ,waveIndex,parentPath) {
  filePath %>%
    getRandomWaveFilePath(parentPath,dataClass, parentPath ) %>%
    readWaveFromNpy(waveIndex)   %>%
    ggplotMFCC( waveIndex, dataClass )
}


plotMFCC3Waves <- function(filePath,dataClass = 1, parentPath = "") {
  plot1 <- filePath %>% plotOneWaveMFCC( dataClass, 1, parentPath )
  plot2 <- filePath %>% plotOneWaveMFCC( dataClass, 2, parentPath )
  plot3 <- filePath %>% plotOneWaveMFCC( dataClass, 3, parentPath )

  ggpubr::ggarrange(plot1,plot2,plot3,common.legend = TRUE)
}


plotOneWaveMFCC2 <- function(filePath, waveIndex , dataClass, ncep = 13,wl = 512, fbtype="htkmel", dcttype="t3"  ) {
  filePath %>%
    readWaveFromNpy(waveIndex)   %>%
    ggplotMFCC( waveIndex, dataClass, filePath, ncep,wl, fbtype, dcttype )
}

showPlotsInGridEqual <- function( plot1 , plot2 , plot3 ) {
  gA <- plot1 %>% ggplot_build() %>% ggplot_gtable()
  gB <- plot2 %>% ggplot_build() %>% ggplot_gtable()
  gC <- plot3 %>% ggplot_build() %>% ggplot_gtable()

  maxWidth = grid::unit.pmax(gA$widths, gB$widths,gC$widths)

  gA$widths <- as.list(maxWidth)
  gB$widths <- as.list(maxWidth)
  gC$widths <- as.list(maxWidth)

  layo <- rbind(c(1,1,1),
                c(2,2,2),
                c(3,3,3))

  grid.newpage()
  grid.arrange(gA, gB,gC, layout_matrix = layo)
}


plotMFCC3Waves2 <- function(filePath, dataClass, ncep = 13,wl = 512, fbtype="htkmel", dcttype="t3") {

  plot1 <- filePath %>% plotOneWaveMFCC2(  1 , dataClass , ncep ,wl , fbtype , dcttype )
  plot2 <- filePath %>% plotOneWaveMFCC2(  2 , dataClass , ncep ,wl , fbtype , dcttype )
  plot3 <- filePath %>% plotOneWaveMFCC2(  3 , dataClass , ncep ,wl , fbtype , dcttype )

  # ggpubr::ggarrange(plot1,plot2,plot3,common.legend = TRUE )
  showPlotsInGridEqual( plot1, plot2 , plot3)
}

savePicturesAsAnimation <- function(pictures,waveClass) {
  pictures %>%
  lapply( magick::image_read) %>%
    magick::image_join() %>%
    magick:: image_animate( fps = 2) %>%
    magick::image_write(  path = sprintf("otherPictures/animations/animationCoeficientsClass%s.gif", waveClass ))
}

animateMFCCCoefPlots <- function(imagesPath, waveClass = 1) {
  imagesPath %>%
    list.files(full.names = TRUE) %>%
    .[15:20] %>%
    savePicturesAsAnimation(waveClass)
}

calCulateMFCCMatrix <- function(wave ,ncep = 13,wl = 512, fbtype="htkmel" , dcttype="t3") {
  wave %>%
    mfccCoefs2(ncep,wl, fbtype, dcttype) %>%
    .$cepstra
}

getSamleOf <- function(nrSamples = 10, target = 1) {
  "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves() %>%
    filter( target ==  target ) %>%
    sample_n( size = nrSamples ) %>%
    select(-X) %>%
    select(filePath) %>%
    unlist() %>%
    as.character() %>%
    lapply(
      function(path) {
        vec <- path %>% readWaveFromNpy(1) %>%
          calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
          as.vector()
        vec %>% set_names( seq(vec %>% length() ) )


      } ) %>%
    bind_rows() %>%
    prcomp( scale. = T) %>%
    .$x %>%
    .[,1:2]%>% as.tibble() %>% add_column(target = target)


}

getSamleOfMixClasses <- function(nrSamples = 10, waveIndex = 1, maxDim = 2 ) {
  wavesClass1 <-   "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves() %>%
    filter( target ==  1 ) %>%
    sample_n( size = nrSamples ) %>%
    select(-X) %>%
    select(filePath) %>%
    unlist() %>%
    as.character() %>%
    lapply(
      function(path,waveIndex) {
        vec <- path %>% readWaveFromNpy(waveIndex) %>%
          calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
          as.vector()
        vec %>% set_names( seq(vec %>% length() ) )


      },waveIndex ) %>%
    bind_rows()

  wavesClass0 <-   "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves() %>%
    filter( target ==  0 ) %>%
    sample_n( size = nrSamples ) %>%
    select(-X) %>%
    select(filePath) %>%
    unlist() %>%
    as.character() %>%
    lapply(
      function(path, waveIndex ) {
        vec <- path %>% readWaveFromNpy(waveIndex ) %>%
          calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
          as.vector()
        vec %>% set_names( seq(vec %>% length() ) )


      }, waveIndex ) %>%
    bind_rows()

  wavesClass0 %>%
    rbind(wavesClass1) %>%
    prcomp( scale. = T) %>%
    .$x %>%
    .[,1:maxDim] %>%
    as.tibble() %>%
    mutate( target = ifelse( row_number() <= nrSamples , 0  , 1  ) %>% as.factor() )



}


getSamleOfMixClassesTsne <- function(nrSamples = 10, waveIndex = 1, maxDim = 2 ) {
  wavesClass1 <-   "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves() %>%
    filter( target ==  1 ) %>%
    sample_n( size = nrSamples ) %>%
    select(-X) %>%
    select(filePath) %>%
    unlist() %>%
    as.character() %>%
    lapply(
      function(path,waveIndex) {
        vec <- path %>% readWaveFromNpy(waveIndex) %>%
          calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
          as.vector()
        vec %>% set_names( seq(vec %>% length() ) )


      },waveIndex ) %>%
    bind_rows()

  wavesClass0 <-   "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves() %>%
    filter( target ==  0 ) %>%
    sample_n( size = nrSamples ) %>%
    select(-X) %>%
    select(filePath) %>%
    unlist() %>%
    as.character() %>%
    lapply(
      function(path, waveIndex ) {
        vec <- path %>% readWaveFromNpy(waveIndex ) %>%
          calCulateMFCCMatrix( ncep = 19,wl = 64, fbtype="htkmel" , dcttype="t3" ) %>%
          as.vector()
        vec %>% set_names( seq(vec %>% length() ) )
      }, waveIndex ) %>%
    bind_rows()

  wavesClass0 %>%
    rbind(wavesClass1) %>%
    scale() %>%
    Rtsne( dims = maxDim) %>%
    .$Y %>%
    as.tibble() %>%
    mutate( target = ifelse( row_number() <= nrSamples , 0  , 1  ) %>% as.factor() )

}

scatterPlot3DDimReduction <- function(data3D) {
  scatterplot3d(x =data3D$V1, y = data3D$V2, z = data3D$V3, color = as.numeric(data3D$target),pch=20,angle=45 ,
                xlab = 'Tsne Dimention x', ylab = 'Tsne Dimention y' , zlab = 'Tsne Dimention z')
}




readFiles <- function() {
  "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves()  %>%
    select(filePath) %>%
    unlist() %>%
    as.character()
}


readLabels <- function() {
  "../machineLearningData/gravitationalWaves/train/file_labels.csv" %>%
    getListOfWaves()  %>%
    select(target) %>%
    unlist()
}

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



dataGeneratorFromFiles1Wave <- function(files, labels,batchSize = 10, waveIndex = 1 ) {
  fileIndex <- 0

  function() {
    if (fileIndex > length(files)) {
      fileIndex <<- 0
    }

    nextIndex = fileIndex + batchSize
    rangeIndexes <- (fileIndex+1): nextIndex
    labelsToprocess <- labels[ rangeIndexes ] %>% as.numeric()
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

#  model %>% compile( loss = "categorical_crossentropy",   optimizer = opt,   metrics = "accuracy"   )
  model %>% compile(
    optimizer = "rmsprop",
    loss = "categorical_crossentropy",
    metrics = c("accuracy")
  )

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



trainWith1Wave <- function(waveIndex = 1){

  labels <- readLabels()
  files <- readFiles()

  gen <- keras:::as_generator.function(dataGeneratorFromFiles1Wave(files, labels, batchSize = 10, waveIndex = 1))

  model <- buildKerasModel1Wave()
  model %>%
    fit_generator( gen,
                  steps_per_epoch = 56000, epochs = 200)
  model

}




