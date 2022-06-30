library(tidyverse)
library(reticulate)
library(Rcpp)
library(RcppCNPy)
library(testthat)


sourceCpp("Functions.cpp")
source_python( "cqtCalculation2.py")
wavePath = "../machineLearningData/gravitationalWaves/train/f/2/f/f2f0bbf138.npy"

cqtCalc = simpleCQT()



test_that("Waves are same 0", {
  expectedWave = wavePath %>% RcppCNPy::npyLoad() %>% .[1,]
  currentWave = wavePath %>% cqtCalc$loadWave(0) %>% as.numeric()
  expect_equal( currentWave , expectedWave )
})

test_that("Waves are same 1", {
  expectedWave = wavePath %>% RcppCNPy::npyLoad() %>% .[2,]
  currentWave = wavePath %>% cqtCalc$loadWave(1) %>% as.numeric()
  expect_equal( currentWave , expectedWave )
})


test_that("Waves are same 2", {
  expectedWave = wavePath %>% RcppCNPy::npyLoad() %>% .[3,]
  currentWave = wavePath %>% cqtCalc$loadWave(2) %>% as.numeric()
  expect_equal( currentWave , expectedWave )
})


test_that("whiten is the same 0", {
  expect_equal( wavePath %>% cqtCalc$loadWave(0)  %>% cqtCalc$whiten() %>% as.numeric() , wavePath %>% cqtCalc$loadWave(0) %>% whiten() %>% as.numeric() )
})


test_that("whiten is the same 1", {
  expect_equal( wavePath %>% cqtCalc$loadWave(1)  %>% cqtCalc$whiten() %>% as.numeric() , wavePath %>% cqtCalc$loadWave(1) %>% whiten() %>% as.numeric() )
})


test_that("whiten is the same 2", {
  expect_equal( wavePath %>% cqtCalc$loadWave(2)  %>% cqtCalc$whiten() %>% as.numeric() , wavePath %>% cqtCalc$loadWave(2) %>% whiten() %>% as.numeric() )
})




test_that("kernals Im is the same", {
  kernalsExpected =  cqtCalc$calcKernels()
  kernalsCurrent = calcKernels()
  expect_equal( calcMatDistance(kernalsCurrent %>% Im(),kernalsCurrent %>% Im() ) %>% rowSums() %>% sum() , 0 )
})


test_that("Kernals Re is the same", {
  kernalsExpected =  cqtCalc$calcKernels()
  kernalsCurrent = calcKernels()
  expect_equal( calcMatDistance(kernalsCurrent %>% Re(),kernalsCurrent %>% Re() ) %>% rowSums() %>% sum(),0 )
})


test_that("hanningWindow padding is the same", {
  wave = wavePath %>% cqtCalc$loadWave(1)
  expected =  wave %>% length() %>%  cqtCalc$getHannWindow()  %>% as.numeric()
  current = wave %>% length() %>% hanningWindow() %>% .[,1] %>% as.numeric()
  expect_equal( expected ,current  )
})



test_that("Signal whitten is the same", {
  expected =  wavePath %>%  cqtCalc$generateWaveWhiten(1) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>% generateWaveWhiten() %>% Re() %>% as.numeric()
  expect_equal( expected ,current )
})


test_that("Signal max is the same", {
  expected =  wavePath %>%  cqtCalc$getMaxWhitenSignal(1)
  current = wavePath %>% cqtCalc$loadWave(1) %>% getMaxWhitenSignal()
  expect_equal( expected ,current )
})



test_that("Signal whitten normalized is the same", {
  expected =  wavePath %>%  cqtCalc$generateWaveWhitenNormalise(1) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>% generateWaveWhitenNormalise() %>% Re() %>% as.numeric()
  expect_equal( expected ,current )
})


test_that("Signal padding is the same", {
  expected =  wavePath %>% cqtCalc$loadWave(1) %>%  cqtCalc$signalPadding(50) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>%  signalPadding(50) %>% as.numeric()
  expect_equal( expected ,current )
})



test_that("Signal padding is the same", {
  expected =  wavePath %>%  cqtCalc$generateWavePadding(1) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>% padding() %>% Re() %>% as.numeric()
  expect_equal( expected ,current )
})


test_that("Signal whiten padding is the same", {
  expected =  wavePath %>% cqtCalc$loadWave(1) %>% cqtCalc$whiten() %>%  cqtCalc$signalPadding(2) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>% whiten() %>% signalPadding(2)  %>% Re() %>% as.numeric()
  expect_equal( expected ,current )
})


test_that("Signal whiten padding is the same", {
  expected =  wavePath %>% cqtCalc$generateWaveWhitenNormalise(1) %>%  cqtCalc$signalPadding(80) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>% generateWaveWhitenNormalise() %>% signalPadding(80)  %>% Re() %>% as.numeric()
  expect_equal( expected ,current )
})



test_that("Signal whiten padding is the same", {
  expected =  wavePath %>% cqtCalc$generateWaveWhitenNormalisePad(1,80) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>% signalPaddingWhitenNormalise(80)  %>% Re() %>% as.numeric()
  expect_equal( expected ,current )
})




test_that("conv1 is the same", {
  expected = cqtCalc$calcConv((c(1:9)  * 1.056) %>% matrix(3,3),c( (1:200) * 1.056 ) ,3) %>% .[1,,]
  current = calcConv((c(1:9) * 1.056)  %>% matrix(3,3), ( c(1:200) * 1.056) ,3) %>% Re()
  expect_equal( expected ,current )
})


# Not Working
test_that("Signal whiten padding is the same", {
  expected =  wavePath %>% cqtCalc$loadWave(1) %>% cqtCalc$generateWhitenPad(20) %>% as.numeric()
  current = wavePath %>% cqtCalc$loadWave(1) %>%  generateWhitenPad(20) %>% Re() %>% as.numeric()
  expect_equal( expected ,current )
})




# Not Working
test_that("conv1 sqt real is the same", {
  expected = wavePath %>% cqtCalc$loadWave(1) %>% cqtCalc$generateCQTConvWaveReal( 24,10 ) %>% .[1,,]
  current = wavePath %>% cqtCalc$loadWave(1) %>% generateCQTConvWaveReal( 24,10 )
  expect_equal( expected ,current )
})



# Not Working
test_that("conv1 sqt imag is the same", {
  expected = wavePath %>% cqtCalc$loadWave(1) %>% cqtCalc$generateCQTConvWaveImag( 24,10 ) %>% .[1,,]
  current = wavePath %>% cqtCalc$loadWave(1) %>% generateCQTConvWaveImag( 24,10 )
  expect_equal( expected ,current )
})





# Not Working
test_that("Signal padding is the same", {

  expected =  wavePath %>%  cqtCalc$generateCQTConvWave(1) %>% .[1,,]
  current = wavePath %>% cqtCalc$loadWave(1) %>% calcCqt() %>% Re()

  expect_equal( expected ,current )
})







test_that("one plus one is two", {
  expect_equal(1 + 1, 2)
})

test_that("you can skip tests if needed", {
  skip("This test hasn't been written yet")
})

test_that("some tests have warnings", {
  expect_equal(log(-1), NaN)
})

test_that("some more successes just to pad things out", {
  expect_true(TRUE)
  expect_false(FALSE)
})