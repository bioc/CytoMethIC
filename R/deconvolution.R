trim <- function(q, rmin=0.15, rmax=0.85) {
    q[q > rmax] <- 1
    q[q < rmin] <- 0
    middle <- na.omit(which(q >= rmin & q <= rmax))
    q[middle] <- (q[middle]-rmin) / (rmax-rmin)
    q
}

## transform fraction (f) by altering 2 components (nu1 and nu2) by step.size
double.update.f <- function(f, nu1, nu2, step.size) {
    f[nu1] <- f[nu1] + step.size
    f[nu2] <- 1 - sum(f[-nu2]) # avoid numerical drift
    if (f[nu1] > 1 || f[nu2] < 0) { return(NULL); }
    else { return(f); }
}

#' Internal function for fraction optimization
#'
#' @param frac initial fraction
#' @param ref reference
#' @param q query
#' @param errFunc error function
#' @param temp annealing temperature
#' @param maxIter maximum iteration to stop after converge
#' @param delta delta score to reset counter
#' @param step.max maximum step, do not adjust
#' @param verbose output debug info
#' @return a list of fractions and min err
.optimizeFrac <- function(frac, ref, q, errFunc,
    temp=0.5, maxIter=1000, delta=0.0001, step.max=1.0, verbose=FALSE) {

    errcurrent <- errFunc(frac, ref, q);
    errmin <- errcurrent
    frac.min <- frac; niter <- 1
    repeat {
        nu <- sample(length(frac), 2) # pick two cell types, including unknown
        step.size <- runif(1) * step.max
        frac.test <- double.update.f(frac, nu[1], nu[2], step.size);
        if (!is.null(frac.test)) {
            if (verbose) {
                message('errcurrent=', errcurrent, 'frac=',
                    paste(lapply(frac, function(x) sprintf('%1.2f', x)),
                        collapse='-'),
                    ';stepsize=', step.size, ';temp=', temp, ';best=',
                    paste(lapply(frac.min,
                        function(x) sprintf('%1.2f', x)), collapse='-'),
                    ';err=', errmin) }

            errtest <- errFunc(frac.test, ref, q)
            if (errtest < errmin) { # update best
                errmin <- errtest
                frac.min <- frac.test
                if ((errmin-errtest) > errmin*delta) { niter <- 1;}
            } else {
                niter <- niter + 1
                if (niter > maxIter) { break;}
            }
            
            ## rejection sampling
            if (runif(1) < exp(-(errtest-errcurrent)/temp)) {
                errcurrent <- errtest
                frac <- frac.test
            }}}

    list(frac=frac.min, err=errmin)
}

#' Reference-based cell type deconvolution
#'
#' This is a reference-based cell composition estimation. The function takes a
#' reference methylation status matrix (rows for probes and columns for cell
#' types) and a query beta value measurement.
#'
#' The length of the target beta values should be the same as
#' the number of rows of the reference Matrix. The function outputs a list
#' containing the estimated cell fraction, the error of optimization.
#' 
#' @param ref reference methylation
#' @param q target measurement: length(q) == nrow(ref)
#' @param trim to trim query input beta values.
#' this relieves unclean background subtraction
#' @param ... extra parameters for optimization.
#' @return a list of fraction, min error.
#' @examples
#'
#' ref = cbind(
#'   CD4 = c(1,1,1,0,1,0),
#'   CD19 = c(0,0,1,1,0,1),
#'   CD14 = c(1,1,1,1,0,1))
#' rownames(ref) = paste0("cg",1:6)
#' trueFrac = runif(3)
#' trueFrac = trueFrac / sum(trueFrac)
#' q = ref %*% trueFrac
#' trueFrac
#' cmi_deconvolution(ref, q)
#' 
#' @export
cmi_deconvolution <- function(
    ref, q, trim=FALSE, ...) {

    if (trim) { q <- trim(q); }
    errFunc <- function(f, ref, q) {
        sum(abs(q - ref %*% f), na.rm=TRUE)
    }

    frac <- rep(1/ncol(ref), ncol(ref))
    res <- .optimizeFrac(frac, ref, q, errFunc, step.max = 0.5, ...)
    res <- .optimizeFrac(res$frac, ref, q, errFunc, step.max = 0.05, ...)

    list(frac = setNames(res$frac, colnames(ref)), err = res$err)
}

#' Reference-based cell type deconvolution (allowing one unknown component)
#'
#' This is a reference-based cell composition estimation. The function takes a
#' reference methylation status matrix (rows for probes and columns for cell
#' types) and a query beta value measurement.
#'
#' The length of the target beta values should be the same as
#' the number of rows of the reference Matrix. The method assumes one unknown
#' component. It outputs a list containing the estimated cell fraction, the
#' error of optimization and methylation status of the unknown component.
#'
#' @param ref reference methylation
#' @param q target measurement: length(q) == nrow(ref)
#' @param trim to trim query input beta values.
#' this relieves unclean background subtraction
#' @param ... extra parameters to .optimizeFrac
#' @return a list of fraction, min error and unknown component methylation state
#' @examples
#'
#' ref = cbind(
#'   CD4 = c(1,1,1,0,1,0),
#'   CD19 = c(0,0,1,1,0,1),
#'   CD14 = c(1,1,1,1,0,1))
#' rownames(ref) = paste0("cg",1:6)
#' trueFrac = runif(4)
#' trueFrac = trueFrac / sum(trueFrac)
#' ref_unk = sample(c(0,1), nrow(ref), replace=TRUE)
#' q = cbind(ref_unk, ref) %*% trueFrac
#' trueFrac
#' res = cmi_deconvolution2(ref, q)
#' res$frac
#' 
#' @export
cmi_deconvolution2 <- function(
    ref, q, trim=FALSE, ...) {

    if (trim) { q <- trim(q); }
    errFunc <- function(f, ref, q) {
        gamma <- q - ref %*% f[2:length(f)]
        sum(ifelse(gamma < f[1] / 2, abs(gamma), abs(gamma - f[1])), na.rm=TRUE)
    }

    frac <- c(1, rep(0, ncol(ref)))
    res <- .optimizeFrac(frac, ref, q, errFunc, step.max=0.5, ...)
    res <- .optimizeFrac(res$frac, ref, q, errFunc, step.max=0.05, ...)
    
    frac <- res$frac
    gamma <- q - ref %*% frac[2:length(frac)]
    g0 <- ifelse(gamma < frac[1] / 2, 0, 1)
    
    list(
        frac = setNames(frac, c("unknown", colnames(ref))),
        err = res$err, g0 = g0)
}
