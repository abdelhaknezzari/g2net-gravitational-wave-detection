library(tidyverse)
11:200 %>% lapply( function( index) {

  nameOfCsv = sprintf( "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.csv",index)
  nameOfRds = sprintf( "../machineLearningData/gravitationalWaves/mfcc/csvMfcc/mfccCoefficients_%s.rds",index)

  nameOfCsv %>%  read.csv() %>% saveRDS(nameOfRds)

} )