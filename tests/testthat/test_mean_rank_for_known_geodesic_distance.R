context("mean_rank_for_known_geodesic_distance")

get_5x5_known_matrix <- function(){
  n_gene = 5
  known = matrix(
    c(0, 1, 1, 0, 0,
      1, 0, 0, 0, 0,
      1, 0, 0, 1, 1,
      0, 0, 1, 0, 0,
      0, 0, 1, 0, 0),
    nrow = n_gene,
    ncol = n_gene,
    byrow = T,
    dimnames = list(sprintf("Gene%s", seq_len(n_gene)),
                    sprintf("Gene%s", seq_len(n_gene))))
  return(known)
}

get_5x5_net_matrix <- function() {
  n_gene = 5
  known = matrix(
    c(0, 5, 3, 4, 0,
      5, 0, 3, 0, 0,
      3, 3, 0, 2, 1,
      4, 0, 2, 0, 2,
      0, 0, 1, 2, 0),
    nrow = n_gene,
    ncol = n_gene,
    byrow = T,
    dimnames = list(sprintf("Gene%s", seq_len(n_gene)),
                    sprintf("Gene%s", seq_len(n_gene)))
  )
  return(known)
}

test_that("mean rank computation", {
  known = get_5x5_known_matrix()
  net = get_5x5_net_matrix()
  ranks = spice::mean_rank_for_known_geodesic_distance(net = net,
                                                       known = known,
                                                       d = 1:2)
  expect_true(ranks[1] == 4.25 && ranks[2] == 5)
})
