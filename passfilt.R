pass.filt1 <- function (y, W, type = c("low", "high", "stop", 
                         "pass"), method = c("Butterworth", "ChebyshevI"), 
          n = 4, Rp = 1) 
{
  if (any(is.na(y))) 
    stop("y contains NA")
  type2 <- match.arg(type)
  nW <- length(W)
  if (type2 == "low" && nW != 1) 
    stop("length(W) > 1")
  if (type2 == "high" && nW != 1) 
    stop("length(W) > 1")
  if (type2 == "stop" && nW != 2) 
    stop("length(W) != 2")
  if (type2 == "pass" && nW != 2) 
    stop("length(W) != 2")
  if (any(W > 1)) {
    f <- 1/W
    p <- W
  }
  else {
    p <- 1/W
    f <- W
  }
  f <- sort(f)
  method2 <- match.arg(method)
  if (method2 == "ChebyshevI") {
    filt <- signal::cheby1(n = n, W = f * 2, type = type2, 
                           Rp = Rp, plane = "z")
  }
  else {
    filt <- signal::butter(n = n, W = f * 2, type = type2, 
                           plane = "z")
  }
  
  yAvg <- mean(y)
  y <- y 
  pad <- max(p) * 2
  ny <- length(y)
  yPad <- c(y[pad:1], y, y[ny:(ny - pad)])
  yFilt <- signal::filtfilt(filt, yPad)
  yFilt <- yFilt[(pad + 1):(ny + pad)]
  return(yFilt) 
}