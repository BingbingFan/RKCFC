# The objective the K-centers iteration targets, and the instrumentation that
# reports how it actually behaves.
#
#     J(c) = sum_i sum_j { Y_ij - Yhat_ij,(c_i) }^2
#
# WHY THIS IS INSTRUMENTED RATHER THAN PROVED
#
# Step 3 minimises J over the labels with the class-specific fits held fixed, so
# the reassignment step is monotone by construction: each subject moves only to
# a class with a strictly smaller residual sum. Step 2 re-estimates each class's
# FPCA under the new partition, and there monotonicity fails. For ordinary
# K-means the analogous step minimises the same objective in closed form, which
# is why Lloyd's algorithm descends. Here the Kendall estimator is not the
# minimiser of a sum of squared reconstruction errors, and local linear
# smoothing of the mean and covariance surfaces carries a bandwidth-dependent
# bias that does not vanish in J. So J can rise across a re-estimation step.
#
# That cannot be argued away, and the same is true of the original KCFC once
# conventional FPCA is replaced. What can be supplied is a measurement of how
# often it happens and by how much, which is what the trace below is for.
#
# THE TRACE IS FREE. univar_clust_update() already returns ss, the residual sum
# of squares of each subject under every class, and it is computed inside the
# reassignment step whether or not anyone looks at it. Summing it two ways gives
# both halves of the decomposition at no additional cost:
#
#   J_before  sum_i ss_i(current label)    after the refit, before reassignment
#   J_after   sum_i min_c ss_i(c)          after reassignment, same fits
#
# J_after <= J_before always, within an iteration; that is Step 3's monotonicity
# and .rkcfc_objective_trace() checks it rather than assuming it. The quantity
# that can rise is J_before at iteration b+1 against J_after at iteration b:
# that difference is the re-estimation step, isolated.

# Plain objective at a given partition, computed from one set of class fits with
# no leave-one-out refitting.
#
# Deliberately independent of the LOO setting, so that two fits are comparable
# whatever produced them. Under LOO = TRUE the ss values inside the loop use a
# leave-one-out fit for a subject's own class, which is the right quantity for
# the assignment RULE but is not sum (Y - Yhat)^2 under a single fit, and is
# therefore not comparable across restarts.
.rkcfc_objective <- function(Ly, Lt, labels, FPCA_C, robust, optns) {
  total <- 0
  for (i in seq_along(Ly)) {
    fit <- FPCA_C[[labels[i]]]
    yhat <- univar_predicted_yhat(fit, Ly[[i]], Lt[[i]], robust, optns)$fittedX
    total <- total + sum((Ly[[i]] - yhat)^2)
  }
  total
}

# Assemble one row of the trace from the reassignment step's own output.
.rkcfc_objective_row <- function(clust_it, current, candidate) {
  ss <- lapply(clust_it, function(x) unlist(x$ss))
  before <- sum(vapply(seq_along(ss), function(i) ss[[i]][current[i]], numeric(1)))
  after  <- sum(vapply(seq_along(ss), function(i) ss[[i]][candidate[i]], numeric(1)))
  c(J_before = before, J_after = after)
}
