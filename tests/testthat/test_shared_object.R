context("SharedObject race condition")
library(parallel)

test_that("sharing objects at the same time", {
  n_gene = 10
  n_sample = 100
  expr = matrix(rnorm(n_gene * n_sample),
                nrow = n_gene,
                ncol = n_sample,
                dimnames = list(sprintf("Gene%s", seq_len(n_gene)),
                                sprintf("Sample%s", seq_len(n_sample))))
  cl <- makeCluster(10)
  tmp <- parallel::clusterEvalQ(cl, {
    # library('spice')
    library(devtools)
    devtools::load_all()
  })
  on.exit(stopCluster(cl))
  doParallel::registerDoParallel(cl)
  expect_error({
    tmp <- foreach::foreach(it = seq_len(50)) %dopar% {
      spice_net = spice::spice(expr, iter = 2, verbose = F, n.cores = 1, seed = 101)
    }
    tmp = gc(verbose = F)
  },NA)
})

