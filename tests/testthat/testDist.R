
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




# Numerical constant
eps <- 10^(-12)

test_that("pqtPareto", {
  
  shape <- 2
  scale <- 3
  
  endpoint <- qpareto(0.9, shape=shape, scale=scale)
  
  ##
  p <- seq(0,1,0.01)
  
  quant <- qtpareto(p, shape=shape, scale=scale, endpoint)
  prob <- ptpareto(quant, shape=shape, scale=scale, endpoint=endpoint)
  
  expect_true( max(abs(prob-p)) < eps )
  
  
  ##
  x <- seq(scale,endpoint,0.01)
  
  prob <- ptpareto(x, shape=shape, scale=scale, endpoint=endpoint)
  quant <- qtpareto(prob, shape=shape, scale=scale, endpoint=endpoint)
  
  expect_true( max(abs(x-quant)) < eps )
  
})



test_that("pqtGPD", {
  
  gamma <- 1/2
  sigma <- sqrt(2)
  mu <- 1.5
  
  
  endpoint <- qgpd(0.9, gamma=gamma, sigma=sigma, mu=mu)
  
  ##
  p <- seq(0,1,0.01)
  
  quant <- qtgpd(p, gamma=gamma, sigma=sigma, mu=mu, endpoint=endpoint)
  prob <- ptgpd(quant, gamma=gamma, sigma=sigma, mu=mu, endpoint=endpoint)
  
  expect_true( max(abs(prob-p)) < eps )
  
  
  ##
  x <- seq(mu,endpoint,0.01)
  
  prob <- ptgpd(x, gamma=gamma, sigma=sigma, mu=mu, endpoint=endpoint)
  quant <- qtgpd(prob, gamma=gamma, sigma=sigma, mu=mu, endpoint=endpoint)
  
  expect_true( max(abs(x-quant)) < eps )
  
})


test_that("pqtBurr", {
  
  alpha <- 2
  rho <- -1
  eta <- 0.5
  
  endpoint <- qburr(0.9, alpha=alpha, rho=rho, eta=eta)
  
  ##
  p <- seq(0,1,0.01)

  quant <- qtburr(p, alpha=alpha, rho=rho, eta=eta, endpoint=endpoint)
  prob <- ptburr(quant, alpha=alpha, rho=rho, eta=eta, endpoint=endpoint)
  
  expect_true( max(abs(prob-p)) < eps )
  
  
  ##
  x <- seq(0,endpoint,0.01)
  
  prob <- ptburr(x, alpha=alpha, rho=rho, eta=eta, endpoint=endpoint)
  quant <- qtburr(prob, alpha=alpha, rho=rho, eta=eta, endpoint=endpoint)
  
  expect_true( max(abs(x-quant)) < eps )

})

test_that("pqtLognormal", {
  
  meanlog <- 2
  sdlog <- 1.5
  
  endpoint <- qlnorm(0.9, meanlog=meanlog, sdlog=sdlog)
  
  ##
  p <- seq(0,1,0.01)
  
  quant <- qtlnorm(p, meanlog=meanlog, sdlog=sdlog, endpoint=endpoint)
  prob <- ptlnorm(quant, meanlog=meanlog, sdlog=sdlog, endpoint=endpoint)
  
  expect_true( max(abs(prob-p)) < eps )
  
  
  ##
  x <- seq(0.01,endpoint,0.01)
  
  prob <- ptlnorm(x, meanlog=meanlog, sdlog=sdlog, endpoint=endpoint)
  quant <- qtlnorm(prob, meanlog=meanlog, sdlog=sdlog, endpoint=endpoint)
  
  expect_true( max(abs(x-quant)) < eps )
  
})


test_that("pqtWeibull", {
  
  shape <- 2
  scale <- 3
  
  endpoint <- qweibull(0.9, shape=shape, scale=scale)
  
  ##
  p <- seq(0,1,0.01)
  
  quant <- qtweibull(p, shape=shape, scale=scale, endpoint)
  prob <- ptweibull(quant, shape=shape, scale=scale, endpoint=endpoint)
  
  expect_true( max(abs(prob-p)) < eps )
  
  
  ##
  x <- seq(scale,endpoint,0.01)
  
  prob <- ptweibull(x, shape=shape, scale=scale, endpoint=endpoint)
  quant <- qtweibull(prob, shape=shape, scale=scale, endpoint=endpoint)
  
  expect_true( max(abs(x-quant)) < eps )
  
})

test_that("pqtExp", {
  
  rate <- pi
  
  endpoint <- qexp(0.9, rate=rate)
  
  ##
  p <- seq(0,1,0.01)
  
  quant <- qtexp(p, rate=rate, endpoint=endpoint)
  prob <- ptexp(quant, rate=rate, endpoint=endpoint)
  
  expect_true( max(abs(prob-p)) < eps )
  
  
  ##
  x <- seq(0,endpoint,0.01)
  
  prob <- ptexp(x, rate=rate, endpoint=endpoint)
  quant <- qtexp(prob, rate=rate, endpoint=endpoint)
  
  expect_true( max(abs(x-quant)) < eps )
  
})