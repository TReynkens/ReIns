
context("Distributions")


test_that("pqSplice", {
  
  # Create MEfit object
  mefit <- MEfit(p=c(0.65,0.35), shape=c(39,58), theta=16.19, M=2)
  
  # Create EVTfit object
  evtfit <- EVTfit(gamma=c(0.76,0.64))
  
  # Create SpliceFit object
  splicefit <- SpliceFit(const=c(0.5,0.996), trunclower=0, t=c(1020,39096), type=c("ME","Pa","Pa"),
                         MEfit=mefit, EVTfit=evtfit)
  
  
  # pSplice(qSplice)
  p <- seq(0,1,0.001)
  quant <- qSplice(p, splicefit)
  prob <- pSplice(quant, splicefit)
  expect_true( max(abs(prob-p),na.rm=TRUE)<10^(-7))
  
  
  # qSplice(pSplice)
  x <- seq(10^2,10^5,10^2)
  prob <- pSplice(x, splicefit)
  quant <- qSplice(prob, splicefit)
  expect_true( max(abs(quant-x)/x,na.rm=TRUE)<10^(-6))

})