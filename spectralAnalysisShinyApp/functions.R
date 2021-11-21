library(RcppCNPy)
library(tidyverse)
library(seewave)
library(tuneR)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)

readWaveFromNpy <- function(filePath,index) {
  filePath %>%
    npyLoad() %>%
    .[index,] %>%
    scale()  %>%
    Wave( left = .,samp.rate = 2048 )

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
