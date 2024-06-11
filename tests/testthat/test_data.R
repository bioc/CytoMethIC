context("data")

test_that("test cmi_predict returns list for RFC", {
  library(sesame)
  library(ExperimentHub)
  library(CytoMethIC)
  betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
  eh <- ExperimentHub()
  modelrfc <- eh[["EH8395"]]
  result <- cmi_predict(imputeBetas(betas, "HM450"), modelrfc)
  expect_is(result, "list")
})


test_that("test cmi_predict returns list for SVM", {
  library(sesame)
  library(ExperimentHub)
  library(CytoMethIC)
  betas <- sesameDataGet("HM450.1.TCGA.PAAD")$betas
  eh <- ExperimentHub()
  modelrfc <- eh[["EH8396"]]
  result <- cmi_predict(imputeBetas(betas, "HM450"), modelrfc)
  expect_is(result, "list")
})

