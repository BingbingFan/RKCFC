# RKCFC 0.2.0

This is the release that produced the results of the accompanying manuscript.
Do not compare results across the 0.1.x boundary, and do not reproduce published
numbers with 0.1.x.

## Bug fixes

* `KCFC_init()` called the non-robust PACE estimator unconditionally, so
  `RKCFC(robust = TRUE)` was started from a non-robust estimator. Since the
  iteration converges to a local optimum, the start selects the basin, so
  results can differ from 0.1.x.
* `silhouette_score()` omitted the covariance bandwidth in its robust branch, so
  the class mean functions behind the coefficient were smoothed differently from
  the model being scored. `robust = FALSE` was unaffected.

## New features

* `RKCFC_restarts()` fits from several random starts and keeps the partition
  with the smallest objective.
* `RKCFC()` gains `init.labels`, to start from user-supplied class labels.
* `KFPCA()` is exported. `scatter = "sign"` gives the spatial-sign covariance in
  place of Kendall's tau, and `meanrobust = TRUE` estimates the mean function by
  a local linear M-estimator with Huber weights. The defaults reproduce the
  previous estimator exactly.
* `normal_scores()` replaces sparse longitudinal measurements by their normal
  scores, ranked among all observations pooled across subjects and times.
* Metrics, preprocessing, and simulation helpers used by the accompanying
  manuscript: `rre()`, `robust_normalize()`, `back_transform()`, `gen_scores()`,
  `vm_intermediate()`, `select_K()`, `ch_index()`, `silhouette_index()`,
  `whiten_scores()`, `spatial_median()`, `sscm()`, `robust_scatter()`, and
  `principal_angles()`.

## Other changes

* The multivariate extension is no longer part of this package and is
  distributed separately as `MRKCFC`.

# RKCFC 0.1.0

* Initial release.
