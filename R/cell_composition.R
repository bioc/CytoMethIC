errFunc <- function(f, g, q) {
    gamma <- q - g %*% f[2:length(f)]
    sum(ifelse(gamma < f[1] / 2, abs(gamma), abs(gamma - f[1])), na.rm=TRUE)
}

## transform fraction (f) by altering 2 components (nu1 and nu2) by step.size
double.transform.f <- function(f, nu1, nu2, step.size) {
    if (f[nu1] + step.size > 1) return(NULL);
    if (f[nu2] - step.size < 0) return(NULL);
    f[nu1] <- f[nu1] + step.size
    f[nu2] <- f[nu2] - step.size
    f[1] <- 1 - sum(f[2:length(f)]) # renormalize to avoid numerical drift
    f
}

dichotomize <- function(q, rmin=0.15, rmax=0.85) {
    q[q > rmax] <- 1
    q[q < rmin] <- 0
    middle <- na.omit(which(q >= rmin & q <= rmax))
    q[middle] <- (q[middle]-rmin) / (rmax-rmin)
    q
}

getg0 <- function(f, g, q) {
    gamma <- q - g %*% f[2:length(f)]
    ifelse(gamma < f[1] / 2, 0, 1)
}

## @param frac0 initial fraction
.optimizeCellComposition <- function(
    g, q, frac0=NULL, temp=0.5, maxIter=1000,
    delta=0.0001, step.max=1.0, verbose=FALSE) {

    M <- ncol(g) # number of reference
    if (is.null(frac0)) { frac <- c(1, rep(0, M))
    } else { frac <- frac0 } # use given fraction estimate
    
    ## initialize
    errcurrent <- errFunc(frac, g, q); errmin <- errcurrent
    frac.min <- frac; niter <- 1
    repeat {
        nu <- sample(seq_len(M+1), 2)
        step.size <- runif(1) * step.max
        frac.test <- double.transform.f(frac, nu[1], nu[2], step.size);
        if (!is.null(frac.test)) {
            if (verbose) {
                message('errcurrent=', errcurrent, 'frac=',
                    paste(lapply(frac, function(x) sprintf('%1.2f', x)),
                        collapse='-'),
                    ';stepsize=', step.size, ';temp=', temp, ';best=',
                    paste(lapply(frac.min,
                        function(x) sprintf('%1.2f', x)), collapse='-'),
                    ';err=', errmin) }

            errtest <- errFunc(frac.test, g, q)
            if (errtest < errmin) { # update best
                errmin <- errtest
                frac.min <- frac.test
                if ((errmin-errtest) > errmin*delta)
                    niter <- 1;
            } else {
                niter <- niter + 1
                if (niter > maxIter) {
                    break;
                }}
            
            ## rejection sampling
            if (runif(1) < exp(-(errtest-errcurrent)/temp)) {
                errcurrent <- errtest
                frac <- frac.test
            }}}

    list(frac.min=frac.min, errmin=errmin)
}

#' Estimate cell composition using reference
#'
#' This is a reference-based cell composition estimation. The function takes a
#' reference methylation status matrix (rows for probes and columns for cell
#' types, can be obtained by getRefSet function) and a query beta value
#' measurement. The length of the target beta values should be the same as
#' the number of rows of the reference matrix. The method assumes one unknown
#' component. It outputs a list containing the estimated cell fraction, the
#' error of optimization and methylation status of the unknown component.
#'
#' @param g reference methylation
#' @param q target measurement: length(q) == nrow(g)
#' @param dichotomize to dichotomize query beta value before estimate,
#' this relieves unclean background subtraction
#' @param refine to refine estimate, takes longer
#' @param ... extra parameters for optimization, this includes
#' temp - annealing temperature (0.5)
#' maxIter - maximum iteration to stop after converge (1000)
#' delta - delta score to reset counter (0.0001)
#' verbose - output debug info (FALSE)
#' @return a list of fraction, min error and unknown component methylation state
## @examples
## g <- diffRefSet(getRefSet(platform='HM450'))
## M <- ncol(g)
## trueFrac <- runif(M+1)
## trueFrac <- trueFrac / sum(trueFrac)
## g0 <- sample(c(0,1), nrow(g), replace=TRUE)
## q <- cbind(g0, g) %*% trueFrac + rnorm(length(g0), mean=0, sd = 0.0)
## q[q<0] <- 0
## q[q>1] <- 1
## est <- estimateCellComposition(g, q)
## @export - TODO fix colinearity
estimateCellComposition <- function(
    g, q, refine=TRUE, dichotomize=FALSE, ...) {

    if (dichotomize) {
        q <- dichotomize(q);
    }

    ## raw
    res <- .optimizeCellComposition(g, q, step.max=1, ...);
    ## refine
    res <- .optimizeCellComposition(
        g, q, frac0=res$frac.min, step.max=0.05, ...)

    list(
        frac = setNames(res$frac.min, c("unknown", colnames(g))),
        err = res$errmin,
        g0 = getg0(res$frac.min, g, q))
}
