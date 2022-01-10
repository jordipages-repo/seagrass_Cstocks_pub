mcheck <- function (obj, ... ) {  
  rs<-resid(obj)
  fv<-fitted(obj)
  par(mfrow=c(1,2))
  require(car)
  plot(fv,rs,xlab="Fitted values",ylab="Residuals")
  abline(h=0, lty=2)
  lines(smooth.spline(fv, rs), col = "red")
  qqPlot(rs,xlab="Normal scores",ylab="Ordered residuals")
  par(mfrow=c(1,1))
  invisible(NULL)}

mcheck2 <- function (obj, ... ) {  
  rs <- resid(obj)
  rs2 <- resid(obj, type = "pearson")
  fv<-fitted(obj)
  par(mfrow=c(2,2))
  require(car)
  plot(fv,rs,xlab="Fitted values",ylab="Residuals")
  abline(h=0, lty=2)
  qqPlot(rs,xlab="Normal scores",ylab="Ordered residuals")
  plot(fv,rs2,xlab="Fitted values",ylab="Residuals Pearson")
  par(mfrow=c(1,1))
  invisible(NULL)}