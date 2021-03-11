context("rank.types")

test_that("spice with rank.types='C' should produces same rank as correlation", {
  n_gene = 10
  n_sample = 100
  expr = matrix(rnorm(n_gene * n_sample),
                nrow = n_gene,
                ncol = n_sample,
                dimnames = list(sprintf("Gene%s", seq_len(n_gene)),
                                sprintf("Sample%s", seq_len(n_sample))))
  spice_net = spice(
    expr,
    iter = 1,
    verbose = F,
    weight.method = "inverse.rank",
    seed = 101,
    frac.gene = 1,
    frac.sample = 1,
    rank.types = "C"
  )
  wgcna_net = abs(cor(t(expr)))
  rs = rank(spice_net[lower.tri(spice_net)])
  rw = rank(wgcna_net[lower.tri(wgcna_net)])
  expect_equal(rs, rw)
  
})
