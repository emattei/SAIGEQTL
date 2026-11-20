# Runs step 1 and step 2 on the bundled regression data and checks the numbers.
# Reference values come from this code at the commit that added the test;
# without a GRM they do not depend on the variance-ratio markers.
# isStoreSigma = TRUE as in the step 1 script; step 2 fails on a model fitted without it.

input <- file.path("..", "..", "extdata", "input")

test_that("step 1 and step 2 reproduce the regression-data results", {
  skip_on_cran()
  skip_if(!dir.exists(input), "extdata/input not available")
  out <- file.path(tempdir(), "regfit")
  on.exit(unlink(paste0(out, "*")), add = TRUE)
  set.seed(1)
  invisible(capture.output(fitNULLGLMM_multiV(
    phenoFile = file.path(input, "seed_1_100_nindep_100_ncell_100_lambda_2_tauIntraSample_0.5_Poisson.txt"),
    phenoCol = "gene_1", covarColList = c("X1", "X2", "pf1", "pf2"), sampleCovarCol = c("X1", "X2"),
    sampleIDColinphenoFile = "IND_ID", traitType = "count", outputPrefix = out,
    isRemoveZerosinPheno = FALSE, isCovariateOffset = FALSE, isCovariateTransform = TRUE, tol = 1e-5,
    plinkFile = file.path(input, "n.indep_100_n.cell_1_01.step1"), IsOverwriteVarianceRatioFile = TRUE,
    useSparseGRMtoFitNULL = FALSE, useGRMtoFitNULL = FALSE, isStoreSigma = TRUE
  )))
  e <- new.env()
  load(paste0(out, ".rda"), envir = e)
  m <- e$modglmm
  expect_true(m$converged)
  expect_equal(sum(m$fitted.values), sum(m$y), tolerance = 1e-6)
  expect_equal(m$theta, c(1, 0.3819597363), tolerance = 1e-6)

  bed <- file.path(input, "n.indep_100_n.cell_1")
  region <- paste0(out, ".region.txt")
  writeLines("2\t1\t9810000", region)
  invisible(capture.output(SPAGMMATtest(
    bedFile = paste0(bed, ".bed"), bimFile = paste0(bed, ".bim"), famFile = paste0(bed, ".fam"),
    chrom = "2", SAIGEOutputFile = paste0(out, ".step2.txt"), minMAF = 0, min_MAC = 20, LOCO = FALSE,
    GMMATmodelFile = paste0(out, ".rda"), varianceRatioFile = paste0(out, ".varianceRatio.txt"),
    SPAcutoff = 2, rangestoIncludeFile = region, markers_per_chunk = 10000
  )))
  res <- read.delim(paste0(out, ".step2.txt"))
  expect_equal(nrow(res), 70)
  expect_equal(res$Tstat[1:5], c(-4.65196, -1.88422, 1.43703, 14.2354, 6.65275), tolerance = 1e-5)
})
