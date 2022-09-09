## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(6, 4)
)

## ----setup--------------------------------------------------------------------
library( binfunest)
set.seed( 31394)

## -----------------------------------------------------------------------------
plot( QPSKdB, 3, 20, log="y", ylim=c( 1e-10, 0.5), panel.first = grid(),
      main="Modulation Performance Curves", xlab="SNR in dB",
      ylab="BER")
curve( QAMdB.8.star, 3, 20, col="blue", add=TRUE)
curve( PSQPSKdB, 3, 20, col="red", add=TRUE)
curve( QAMdB.16, 3, 20, col="green", add=TRUE)
legend( "bottomleft", legend=c( "QPSKdB", "QAMdB.8.star", "PSQPSKdB",
                                "QAMdB.16"),
        lty=c( 1, 1, 1, 1),
        col=c("black", "blue", "red", "green"))

## ----curves, fig.cap="Dotted line shows B2B limiting BER."--------------------
QPSKdB.B2B <- B2BConvert( QPSKdB)
O1 <- 3
B1 <- 16
s <- 0:20
plot( s, y=QPSKdB( s), typ="l", log="y", ylim=c( 1e-10, 0.5), 
      panel.first = grid(), main="QPSK Performance Limits", xlab="SNR in dB",
      ylab="BER")
lines( s, y=QPSKdB.B2B( s, Inf, O1), col="blue", lwd=2)
lines( s, y=QPSKdB.B2B( s, B1, O1), col='red')
lines( s, y=QPSKdB.B2B( rep( B1, length( s)), Inf, O1), col='black', lty=2)
legend( "bottomleft", legend=c( "Theory", "Offset", "Offset + B2B", "B2B"),
        lty=c( 1, 1, 1, 2), lwd=c(1, 2, 1, 1), 
        col=c( 'black', 'blue', 'red', 'black'))

## ----data, fig.cap="Samples above 17 dB resulted in zero errors (in this example) and so are not plotted."----
N <- 1000000
(r <- rbinom( length( s), N, QPSKdB.B2B( s, B1, O1)))
bplot1 <- function( s, r, N, O, B, ylim=c( 1e-10, 0.5)) {
  plot( s, y=QPSKdB.B2B( s, B, O), col='red', log='y', type='l', ylim=ylim,
        main="QPSK Performance Limits", xlab="SNR in dB",
      ylab="BER", panel.first = grid())
  points( s, r/N)
  lines( s, y=QPSKdB( s))
  lines( s, y=QPSKdB.B2B( s, Inf, O), col="blue", lwd=2)
  lines( s, y=QPSKdB.B2B( rep( B1, length( s)), Inf, O1), col='black', lty=2)
  legend( "bottomleft",
        legend=c( "Data", "Theory", "Offset", "Offset + B2B", "B2B"),
        lty=c( NA, 1, 1, 1, 2), lwd=c(NA, 1, 2, 1, 1), pch=c( 1, NA, NA, NA, NA),
        col=c( 'black', 'black', 'blue', 'red', 'black'))
}
bplot1( s, r, N, O1, B1)

## ---- fig.cap="Non-linear least squares fit (green line) of the line to the data."----
df1 <- data.frame( Errors=r, Trials=rep( N, length( s)), SNR=s)
nls.fit1 <- nls( Errors/Trials ~ QPSKdB.B2B( SNR, b, o), data=df1, 
                 start=list( b=10, o=0))
summary( nls.fit1)
bplot1(s, r, N, O1, B1)
c <- coef( nls.fit1)
lines( s, y=QPSKdB.B2B( s, c['b'], c['o']), col='green')

## ---- fig.cap="Non-linear fit (green line) excluding the zeros."--------------
nls.fit2 <- nls( Errors/Trials ~  QPSKdB.B2B( SNR, b, o), 
                 data=subset( df1, Errors > 0), start=list( b=10, o=0))
summary( nls.fit2)
bplot1(s, r, N, O1, B1)
c <- coef( nls.fit2)
lines( s, y=QPSKdB.B2B( s, c['b'], c['o']), col='green')

## ---- fig.cap="Linear plot of performance showing how close to zero the estimates lie."----
 plot( s, y=QPSKdB.B2B( s, B1, O1), col='red', type='l', panel.first = grid(),
       main="QPSK Performance Limits", xlab="SNR in dB", ylab="BER")
  points( s, r/N, panel.first = grid())
  lines( s, y=QPSKdB( s))
  lines( s, y=QPSKdB.B2B( s, Inf, O1), col="blue", lwd=2)
  lines( s, y=QPSKdB.B2B( s, c['b'], c['o']), col='green')
  legend( "bottomleft",
        legend=c( "Data", "Theory", "3 dB Offset", "Offset + B2B", "Fit"),
        lty=c( NA, 1, 1, 1, 1), lwd=c(NA, 1, 2, 1, 1), pch=c( 1, NA, NA, NA, NA),
        col=c( 'black', 'black', 'blue', 'red', 'green'))


## ---- fig.cap="Results of non-linear fit (green line) to Q function.  Zeros are still filtered out."----
nls.fit3 <- nls( Q_( Errors/Trials) ~  Q_( QPSKdB.B2B( SNR, b, o)), 
                 data=subset( df1, Errors > 0), start=list( b=10, o=0))
summary( nls.fit3)
bplot1(s, r, N, O1, B1)
c <- coef( nls.fit3)
lines( s, y=QPSKdB.B2B( s, c['b'], c['o']), col='green')

## -----------------------------------------------------------------------------
dbinom( r, N, QPSKdB.B2B( s, B1, O1), log=TRUE)

## ----llsb---------------------------------------------------------------------
llsb <- function( par) 
  -sum( sort( dbinom( r, N, QPSKdB.B2B( s, par[1], par[2]), log=TRUE),
              decreasing=TRUE))
mle1 <- optim( par=c( b2b=20, off=0), llsb, hessian=TRUE)
mle1$par
(mle1s <- sqrt( diag( solve( mle1$hessian))))

## -----------------------------------------------------------------------------
c( (mle1$par["b2b"] - B1)/mle1s["b2b"], (mle1$par["off"] - O1)/mle1s["off"])

## ----mlef---------------------------------------------------------------------
mlef <- function( N, s, r, f, ..., startpar=c( b2b=20, off=0),
                    method="Nelder-Mead", lower, upper, control) {
  stopifnot( length( N) == length( s) && length( s) == length( r))
  ll <- function( par) {
    -sum( dbinom( r, N, f( s, par[1], par[2], ...), log=TRUE))
  }
  res <- optim( startpar, ll, hessian=TRUE)
}

mle2 <- mlef( rep( N, length( s)), s, r, QPSKdB.B2B)
mle2$par
sqrt( diag( solve( mle2$hessian)))

## ---- fig.cap="Estimated as green dashed line (dashing allows the red line to peek through)."----
bplot1(s, r, N, O1, B1)
lines( s, QPSKdB.B2B( s, mle2$par[1], mle2$par[2]), col='green', lty=4)
legend( "bottomleft",
        legend=c( "Data", "Theory", "3 dB Offset", "Offset + 15dB B2B",
                  "Estimated"),
        lty=c( NA, 1, 1, 1, 4), pch=c( 1, NA, NA, NA, NA),
        col=c( 'black', 'black', 'blue', 'red', 'green'))

## ----stats4_mle---------------------------------------------------------------
require( stats4)
llsb2 <- function( b2b, off) 
   -sum( dbinom( r, N, QPSKdB.B2B( s, b2b, off), log=TRUE))
mle3 <- mle( llsb2, start=c( b2b=20, off=0), nobs=length(s))
summary( mle3)

## ----loglik-------------------------------------------------------------------
logLik( mle3)
vcov( mle3)
plot(profile( mle3), absVal = FALSE)
confint( mle3)


## ----mleB2B-------------------------------------------------------------------
df <- data.frame( Errors=r, SNR=s)
(est1 <- mleB2B( data=df, Errors="Errors", N=N, f=QPSKdB.B2B, 
                 fparms=list( x="SNR"), start=c(b2b=20, offset=0)))

## ----vhB2B, fig.cap="24 dB B2B Q and 3 dB Offset estimates."------------------
O2 <- 3
B2 <- 24
N <- 1000000
(r2 <- rbinom( length( s), N, QPSKdB.B2B( s, B2, O2)))
mle4 <- mleB2B(  Errors=r2, N=N, f=QPSKdB.B2B,
                 fparms=list( x=s), start=c(b2b=20, offset=0))
summary(mle4)
mle4coef <- coef(mle4)
plot(  s, r2/N, log='y', panel.first = grid(), ylim=c(1e-14, 0.5))
lines( s, y=QPSKdB( s), col='black')
lines( s, y=QPSKdB.B2B( s, Inf, O1), col="blue", lwd=2)
lines( s, y=QPSKdB.B2B( s, B2, O2), col="red")
lines( s, y=QPSKdB.B2B( s, mle4coef[1],  mle4coef[2]), col="green")
legend( "bottomleft",
        legend=c( "Data", "Theory", "Offset", "Offset + B2B",  
                  "Estimated"),
        lty=c( NA, 1, 1, 1), lwd=c(NA, 1, 2, 1, 1), 
        col=c( 'black', 'black', 'blue', 'red', 'green'),
        pch=c( 1, NA, NA, NA, NA))

## ----mle5---------------------------------------------------------------------
(mle5 <- mleB2B( Errors=r2, N=N, f=QPSKdB.B2B, method="Brent", 
                 fparms=list( x=s, B2B=+Inf), start=c( offset=0), 
                 lower=-6, upper=10))

## -----------------------------------------------------------------------------
AIC(mle4, mle5)

## ----mle6, fig.cap="24 dB B2B Q and 3 dB Offset estimates with additional data point."----
N2 <- 100e9*3600
Nvec <- c( rep(N, length(s)), N2)
s2 <- c(s, 19)
(r4 <- c( r2, rbinom( 1, N2, QPSKdB.B2B( 19, B2, O2))))
mle6 <- mleB2B(  Errors=r4, N=Nvec, f=QPSKdB.B2B,
                 fparms=list( x=s2), start=c(b2b=20, offset=0))
summary( mle6)
mle6coef <- coef(mle6)
plot(  s2, r4/Nvec, log='y', panel.first = grid(), ylim=c(1e-14, 0.5))
lines( s2, y=QPSKdB( s2), col='black')
lines( s2, y=QPSKdB.B2B( s2, Inf, O1), col="blue", lwd=2)
lines( s2, y=QPSKdB.B2B( s2, B2, O2), col="red")
lines( s2, y=QPSKdB.B2B( s2, mle6coef[1],  mle6coef[2]), col="green")
legend( "bottomleft",
        legend=c( "Data", "Theory", "Offset", "Offset + B2B",  
                  "Estimated"),
        lty=c( NA, 1, 1, 1), lwd=c(NA, 1, 2, 1, 1), 
        col=c( 'black', 'black', 'blue', 'red', 'green'),
        pch=c( 1, NA, NA, NA, NA))

## ----noB2B--------------------------------------------------------------------
(r3 <- rbinom( length( s), N, QPSKdB( s - O2)))
mle7 <- mleB2B( N=N, Errors=r3, f=QPSKdB.B2B, fparms=list( x=s),
                start=c(b2b=20, offset=0))
summary(mle7)
mle8 <- mleB2B( N=N, Errors=r3, f=QPSKdB.B2B, fparms=list( x=s, B2B=+Inf),
                start=list(offset=0), method="Brent", lower=-6, upper=10)
summary(mle8)


AIC(mle7,mle8)
exp( (AIC( mle7) - AIC( mle6)) / 2)

## ----Qab----------------------------------------------------------------------
Q_ab <- function( s, a, b) Q_( sqrt( a^2 * s/(1 + b^2 * s)))
mle8 <- mleB2B( N=N, Errors=r2, f=Q_ab, start=c( a=1, b=0), fparms=list(s=undB(s)))
summary( mle8)

## ----plotQab------------------------------------------------------------------
mle8coef <- coef( mle8)
plot( s, y=r2/N, log='y', type='p', panel.first = grid())
lines( s, QPSKdB( s))
lines( s, Q_(sqrt( undB( s)/(1 + undB( -80) * s))), col='red')
lines( s, y=Q_ab( undB( s), mle8coef[1],  mle8coef[2]), col="green")
legend( "bottomleft",
        legend=c( "Data",  "Theory", "0 dB Offset + 80 dB B2B", 
                  "Estimated"),
        lty=c( NA, 1, 1, 1), col=c( 'black', 'black', 'red', 'green'),
        pch=c( 1, NA, NA, NA))

## -----------------------------------------------------------------------------
O4 <- 3
B4 <- 15
(r4 <- rbinom( length( s), N, QPSKdB.B2B( s, B4, O4)))
mle9 <- mleB2B( N=N, Errors=r4, f=Q_ab, start=c( a=1, b=0), 
                fparms=list(s=undB(s)))
summary( mle9)
mle9coef <- coef( mle9)
(mle9sd <- sqrt( diag( vcov( mle9))))
plot( s, r4/N, log='y',panel.first = grid())
lines( s, QPSKdB( s))
lines( s, y=QPSKdB.B2B( s, Inf, O4), col="blue")
lines( s, y=QPSKdB.B2B( s, B4, O4), col="red")
lines( s, y=Q_ab( undB( s), mle9coef[1],  mle9coef[2]), col="green", lty=5)
legend( "bottomleft",
        legend=c( "Data", "Theory", "Theory + 3 dB", "3 dB Offset + 20 dB B2B",
                  "Estimated"),
        lty=c( NA, 1, 1, 1, 5), 
        col=c( 'black', 'black', 'blue', 'red', 'green'),
        pch=c( 1, NA, NA, NA, NA))

## -----------------------------------------------------------------------------
-2*dB(mle9coef)

## -----------------------------------------------------------------------------
Q8Sbits <- BERDFc[BERDFc$Name=="8QAM Star",] # Select one constellation
Q8Sbits$Nbits <- Q8Sbits$Bps * Q8Sbits$N # Add  a Nbits column
QAMdB.8.star.B2B <- B2BConvert( QAMdB.8.star) 
mle10 <- mleB2B( data=Q8Sbits, N="Nbits", Errors="BER", 
                 f=QAMdB.8.star.B2B , start=c( B2B=30, offset=0), 
                 fparms=list( x="SNR"))
summary( mle10)

## -----------------------------------------------------------------------------
mle11 <- mleB2B( data=Q8Sbits, N="Nbits", Errors="BER", 
                 f=QAMdB.8.star.B2B , method="Brent", 
                 fparms=list( x="SNR", B2B=+Inf), start=c( offset=0), 
                 lower=-6, upper=10)
summary( mle11)
AIC( mle10, mle11)

## -----------------------------------------------------------------------------
Q8Sbits

