#' Master data frame for all model objects
#'
#' This is an internal object which will be updated on every new release
#'
#' @name cmi_models
#' @docType data
#' @format tibble
#' @return master sheet of CytoMethIC model objects
#' @examples print(cmi_models[,c("EHID","Title")])
#' @export
NULL

## internal function for checking features
clean_features <- function(
    betas, cmi_model, verbose = FALSE) {
    
    features <- cmi_model$features
    idx <- match(features, names(betas))
    if (sum(is.na(idx)) == 0 && sum(is.na(betas[idx])) == 0) {
        return(betas[idx])
    }
    
    if (verbose) {
        warning(sprintf("Missing %d/%d features.",
            sum(is.na(betas[idx])), length(idx)))
    }

    ## if still have missing values
    idx <- match(features, names(betas))
    if (sum(!is.na(idx)) == 0) {
        stop(sprintf("No overlapping probes.\n%s\n%s\n",
            "If on different platform, use mLiftOver.",
            "If missing value, use imputeBetas"))
    }
    if (sum(!is.na(betas[idx])) == 0 ||  # all-NA
        (sum(is.na(betas[idx])) > 0 && ( # some NA and model requires all
            is.null(cmi_model$features_require_all) ||
            cmi_model$features_require_all))) {
        stop(sprintf("Missing %d/%d features.\n%s\n%s\n",
            sum(is.na(betas[idx])), sum(!is.na(idx)),
            "If on different platform, use mLiftOver.",
            "If missing value, use imputeBetas."))
    }
    betas <- betas[features]
}

#' The cmi_predict function takes in a model and a sample, and uses the model
#' to predict it.  This function supports randomForest, e1071::svm, xgboost,
#' and keras/tensorflow models. For xgboost and keras models, the features used
#' in classification as well as a label mapping must be provided for output.
#' 
#' @param betas DNA methylation beta
#' @param cmi_model Cytomethic model downloaded from ExperimentHub
#' @param BPPARAM use MulticoreParam(n) for parallel processing
#' @param verbose be verbose with warning
#' @return predicted cancer type label
#' @examples
#' 
#' library(sesame)
#' library(ExperimentHub)
#' library(CytoMethIC)
#'
#' ## Cancer Type
#' model = ExperimentHub()[["EH8395"]]
#' betas = openSesame(sesameDataGet("EPICv2.8.SigDF")[[1]])
#' betas = imputeBetas(mLiftOver(betas, "HM450"))
#' cmi_predict(betas, model)
#' 
#' betas = openSesame(sesameDataGet('EPIC.1.SigDF'), mask=FALSE)
#' cmi_predict(betas, model)
#'
#' betas = sesameDataGet("HM450.1.TCGA.PAAD")$betas
#' betas = imputeBetas(betas)
#' cmi_predict(betas, model)
#' 
#' @import stats
#' @import tools
#' @import sesameData
#' @import sesame
#' @import ExperimentHub
#' @import BiocParallel
#' @importFrom methods is
#' @export
cmi_predict <- function(betas, cmi_model,
    verbose = FALSE, BPPARAM = SerialParam()) {
    
    if(names(cmi_model)[[1]] == "model_serialized") {
        requireNamespace("keras")
        requireNamespace("tensorflow")
        cmi_model <- list(
            model = keras::unserialize_model(cmi_model$model_serialized),
            features = cmi_model$features,
            label_levels = cmi_model$label_levels)
    }

    if (is.matrix(betas)) {
        res <- bplapply(seq_len(ncol(betas)), function(i) {
            cmi_predict(betas[,i], cmi_model, verbose = verbose)
        }, BPPARAM = BPPARAM)
        names(res) <- colnames(betas)
        return(res)
    }

    stopifnot(is.numeric(betas))
    betas <- clean_features(betas, cmi_model, verbose = verbose)

    if (is(cmi_model$model, "function")) {
        cmi_model$model(betas, cmi_model$model_data)
    } else if (is(cmi_model$model, "randomForest")) {
        requireNamespace("randomForest")
        res <- sort(predict(cmi_model$model,
            newdata = t(betas), type = "prob")[1,], decreasing = TRUE)
        list(response = names(res)[1], prob = res[1])
    } else if (is(cmi_model$model, "svm")) {
        betas <- t(as.data.frame(betas))
        if (!requireNamespace("e1071", quietly = TRUE)) stop("e1071 not installed")
        model <- cmi_model$model
        betas <- t(as.data.frame(betas[, attr(model$terms, "term.labels")]))
        res <- as.character(predict(model, newdata = betas))
        probs <- attr(predict(model, newdata = betas, probability = TRUE), "probabilities")
        prob_max <- apply(probs, MARGIN = 1, FUN = max)[1]
        list(response = as.character(res), prob = prob_max)
    } else if (is(cmi_model$model, "xgb")) {
        betas <- t(as.data.frame(betas))
        if (!requireNamespace("xgboost", quietly = TRUE)) stop("xgboost not installed")
        if (is.null(cmi_model$features)) stop("Must provide feature parameter with xgboost model")
        if (is.null(cmi_model$label_levels)) stop("Must provide label_levels parameter with xgboost model")
        betas <- (t(as.data.frame(betas[, cmi_model$features])))
        model <- cmi_model$model
        betas <- xgboost::xgb.DMatrix(t(as.matrix(betas)))
        pred_probabilities <- predict(model, betas)
        num_classes <- length(pred_probabilities)
        pred_prob_matrix <- matrix(pred_probabilities, nrow = 1, ncol = num_classes, byrow = TRUE)
        max_probability <- apply(pred_prob_matrix, 1, max)
        pred_label <- cmi_model$label_levels[apply(pred_prob_matrix, 1, which.max)]
        list(response = pred_label, prob = max_probability)
    } else if (is(cmi_model$model, "keras")) {
        if (is.null(cmi_model$features)) stop("Must provide feature parameter with Keras model")
        if (is.null(cmi_model$label_levels)) stop("Must provide label_levels parameter with Keras model")
        betas <- t(as.data.frame(betas))
        betas <- t(as.data.frame(betas[, cmi_model$features]))
        betas <- t(as.matrix(betas))
        pred_prob_matrix <- predict(cmi_model$model, betas)
        max_probability <- apply(pred_prob_matrix, 1, max)
        highest_prob_prediction <- apply(pred_prob_matrix, 1, function(x) which.max(x))
        pred_label <- cmi_model$label_levels[highest_prob_prediction]
        list(response = pred_label, prob = max_probability)
    } else {
        stop("Model not supported by package.")
    }
}



#' Check CytoMethIC versions
#'
#' print package verison of cytomethic and depended packages to help
#' troubleshoot installation issues.
#'
#' @return print the versions of cytomethic and dependencies
#' @importFrom utils packageVersion
#' @importFrom BiocManager version
#' @export
#' @examples
#' cmi_checkVersion()
cmi_checkVersion <- function() {
    rv <- R.Version()
    message(
      "CytoMethIC requires matched versions of ",
      "R, sesame, sesameData and ExperimentHub.\n",
      "Here is the current versions installed:\n",
      "R: ", rv$major, rv$minor, "\n",
      "Bioconductor: ", version(), "\n",
      "CytoMethIC: ", packageVersion("CytoMethIC"), "\n",
      "sesame: ", packageVersion("sesame"), "\n",
      "sesameData: ", packageVersion("sesameData"), "\n",
      "ExperimentHub: ", packageVersion("ExperimentHub"), "\n"
    )
}
