#' Theoretical error rate functions
#'
#' Functions to calculate the theoretical performance of common modulation
#' formats.  Includes the functions `dB (x)` (returns `10log10(x)`), `undB(x)`
#' (reverses `dB(x)`), `Q_( x)` (Markum's Q function), and `Q_Inv(x)`
#' (returns the
#' SNR in Decibels to get probability x).  Also includes `mod_Inv`, which returns
#' the SNR required for a the function `f` to reach the supplied BER (bit
#' error rate, or bit error probability).
#'
#' The rest of the functions return the probability of a bit error given the
#' SNR in Decibels.
#' *  `QPSKdB` is Quadrature Phase shift keyed: two bits per symbol.
#' *  `DQPSK` is differentially detected differentially coded QPSK.
#' *  `DQPSKDDdB` is differentially detected differential QPSK (coherently
#'      detected but differentially decoded. See `DQPSK` above.
#' *  `PSQPSKdB` is polarization-shifted QPSK: it is dual pole, but only
#'      one pole is active at any one time, thus supplying three bits per
#'      symbol. (See Agrell & Karlsson (2009, DOI:10.1109/JLT.2009.2029064)).
#' *  `MPSKdB(x, M)` is generic M-ary phase shift keying of `M` points in a circle.
#' *  `MPSKdB.8` simply returns `MPSKdB(x, 8)`
#' *  `QAMdB.8.star` is the optimal star configuration of 8-ary Quadrature
#'      Amplitude Modulation (QAM), such that
#'      the legs are at \eqn{\pm1} and \eqn{\pm(1+\sqrt3)}.
#' *  `QAMdB(x, M)` is generic rectangular QAM constellation of `M` points.
#' *  `QAMdB.16` Returns the BER for the rectangular QAM constellation according to
#'      Proakis Eq. 5-2-80.
#' * `mod_Inv` will take a function `f(x)` and return the x such that
#'      `f(x)==perr`
#'      but it does this based on the `log( f(x))` and the `log( perr)`, so
#'      `f(x)>0`.
#' * `mod_InvV` is a vectorized version (give it a vector of BERs and it returns a
#'      vector of SNRs).
#'
#' @param x The SNR (Eb/No) in Decibels, possibly a vector.
#' @param M The integer number of symbols > 4.
#' @param perr a probability of a bit error.
#' @param f a function (usually a BER function).
#' @param pv a vector of BERs.
#' @param guess a guess for the `perr` (the default usually works).
#' @param offset an offset in Decibels for guesses in `mod_InvV`.
#'
#' @name Theoretical
NULL
#>

sqrt_2 <- sqrt( 2.0)
tms3_ <- 3.0 - sqrt( 3.0) # Need this for 8QAM Star.

#' @rdname Theoretical
#' @name is.wholenumber
#' @param x a real number
#' @param tol the tolerance to test x with.
#' @return `is.wholenumber(x)` returns `TRUE` if `c-round(x) < tol`.
#' @export
is.wholenumber <- function( x, tol = sqrt( .Machine$double.eps))
   abs(x - round( x)) < tol # Returns TRUE if double x is close to integer.

#' @rdname Theoretical
#' @return \code{dB(x)} returns \code{10*log10(x)}
#' @examples
#' dB( 10) # == 10
#' @export
dB <- function( x) 10.0 * log10( x)

#' @rdname Theoretical
#' @return \code{undB(x)} returns \code{10^(x/10)}
#' @examples
#' undB( 20) # == 100
#' @export
undB <- function( x) 10^( 0.1 * x)

#' @rdname Theoretical
#' @export
Q_ <- function( x) stats::pnorm( x, lower.tail=FALSE)

#' @rdname Theoretical
#' @return \code{Q_Inv(x)} returns \code{2*dB( -qnorm(x))}, which is the
#' SNR (in Decibels) required to get a probability of error of x.
#' Q_Inv( Q_( undB( x/2))) = x and Q_( undB( Q_Inv( x)/2))=x
#'
#' @examples
#' Q_Inv( Q_( undB( 10/2))) # = 10
#' Q_( undB( Q_Inv( 0.001)/2)) # = 0.001
#'
#' @export
Q_Inv <- function( perr) 2.0 * dB( -stats::qnorm( perr))

#' @rdname Theoretical
#' @export
QPSKdB <- function( x) Q_( sqrt_2 * undB( 0.5 * x))

#' @rdname Theoretical
#' @export
DQPSKdB <- function( x) 0.5 * exp( -undB( x))

# Differentially decoded DQPSK.
#' @rdname Theoretical
#' @export
DQPSKDDdB <- function( x) {p <- QPSKdB( x); 2.0 * p * (1.0 - p)}

# Polarization shifted QPSK. Aggrell & Karlson
# (2009, DOI:10.1109/JLT.2009.2029064)
#' @rdname Theoretical
#' @export
PSQPSKdB <- function( x) {
   p <- Q_( sqrt( 3.0 * undB( x))); (7.0 / 6.0) * (2.0 - p) * p}

#' @rdname Theoretical
#' @export
MPSKdB <- function( x, M){
   stopifnot( is.wholenumber( M))
   stopifnot( M > 4) # 4 is QPSK, 2 is BPSK
   k <- log2( M);
   (2.0 / k) * Q_( sqrt( 2.0 * k * undB( x)) * sin( pi / M))
}

#' @rdname Theoretical
#' @export
MPSKdB.8 <- function (x) MPSKdB( x, 8)

#' @rdname Theoretical
#' @export
QAMdB.8.star <- function( x) 1.25 * Q_( sqrt( tms3_ * undB( x)))

#' @rdname Theoretical
#' @export
QAMdB <- function( x, M) { # M-ary QAM
   stopifnot( is.wholenumber( M))
   stopifnot( M > 4)
   k <- log2( M)
   (4 / k) * Q_( sqrt( 3.0 * k * undB(x) / (M - 1)))
}

#' @rdname Theoretical
#' @export
QAMdB.16 <- function( x) (1.0 - (1.0 - 1.5 * Q_( sqrt( 0.8 * undB( x))))^2) / 4
#Proakis Eq. 5-2-80

#' @rdname Theoretical
#' @returns \code{ mod_Inv( f, x)} returns a list with the SNR in Decibels to
#' reach the BER
#' \code{perr} such that \code{f( mod_Inv( f, x)$x) = x}.
#' The returned list has elements
#'  \code{$x} as the SNR and
#'  \code{$fval} as the function value.
#'
#' @examples
#' mod_Inv( QPSKdB, QPSKdB( 7)) # yields 7
#'
#' @seealso [pracma::fzero()]
#' @export
## mod_Inv will take a function f(x) and return the x such that f(x)==perr
## but it does this based on the log(f(x)) and the log(perr), so f(x)>0.
mod_Inv <- function( f, perr, guess=Q_Inv( perr))
   if (perr>0)  {pracma::fzero( function( X) log( f( X)) - log( perr), guess)
   } else list( x=Inf, fval=0)

#' @rdname Theoretical
#' @examples
#' mod_InvV(QPSKdB, QPSKdB(c(6,7)))
#'
#' @export
mod_InvV <- function( f, pv, offset=0.0)
   # vectorized version of mod_Inv. pv can be a vector of BERs.
   # The offset is an offset to the Q_ function (in dB), such that
   # the guesses passed to fzero will be Q_Inv( p) - offset.
   apply( cbind( pv, Q_Inv( pv) - offset), 1,
          function( x) mod_Inv( f, x[ 1], x[ 2])$x)
