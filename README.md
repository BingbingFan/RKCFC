# RKCFC

Robust K-centers functional clustering for sparsely and irregularly observed
longitudinal data that may be non-Gaussian.

Trajectory clustering methods built on functional principal component analysis
assume, at some point, that a sample covariance operator is a reliable estimate.
For medical longitudinal measurements that is often false. Serum tumour markers,
for instance, run orders of magnitude above their median in a minority of
patients, and a covariance surface estimated from squared deviations follows
those few subjects. RKCFC replaces the covariance operator inside the K-centers
iteration with Kendall's tau operator, which has the same eigenspace under an
elliptical process but needs no moment assumption on the principal component
scores.

## Installation

```r
# install.packages("remotes")
remotes::install_github("BingbingFan/RKCFC")

# the release that produced the results of the accompanying manuscript
remotes::install_github("BingbingFan/RKCFC@v0.2.0")
```

## Quick start

```r
library(RKCFC)
data(pbc2)

# Serum bilirubin is sparse, irregularly timed, and strongly right-skewed.
pbc  <- pbc2[!is.na(pbc2$serBilir), c("id", "year", "serBilir")]
keep <- names(which(table(pbc$id) >= 3))
pbc  <- pbc[pbc$id %in% keep, ]

interval <- range(pbc$year)
optns <- list(
  interval = interval,
  nK       = 2,                       # components retained
  bw       = 0.1 * diff(interval),    # subject-level smoothing; no GCV path, must be numeric
  bwK      = "GCV",                   # Kendall surface bandwidth
  nRegGrid = 51,
  basis    = fda::create.bspline.basis(
    interval, nbasis = 13, norder = 4,
    breaks = seq(interval[1], interval[2], length.out = 11))
)

fit <- RKCFC(data = pbc, id = "id", ys = "serBilir", timevar = "year",
             K = 2, robust = TRUE, optns = optns,
             init.method = "two-staged", parallel = FALSE, LOO = FALSE)

table(fit$res_clust$g)
```

`robust = FALSE` gives the K-centers functional clustering (KCFC) model of Chiou
and Li (2007), which is the natural comparator, the same iteration with the PACE
covariance estimator in place of Kendall's tau.

## Choosing the number of classes

`silhouette_score()` is the modified silhouette coefficient for sparse designs.
Fit over a range of `K` and take the maximiser.

Two things to know before reporting the result. It is undefined at `K = 1`, so a
value of 0 there is a convention rather than a number. And because both
distances are measured to class *mean functions*, it rewards classes separated in
level and penalises classes separated in shape or slope, so a method producing the
latter can score lower without being worse. Read a small maximum against what the
criterion returns when the number of classes is known, not in isolation.
`select_K()` supplies score-space alternatives for that reason.

## What is robust, and what is not

Only the covariance surface and the eigenfunctions.

| Component | Estimator | Robust |
|---|---|---|
| Kendall's tau covariance surface, eigenfunctions | Kendall's tau operator | **yes** |
| Mean function | local linear least squares | no |
| Principal component scores | least squares | no |
| Subject-level pre-smoothing | Nadaraya–Watson | no |
| Class-assignment rule | squared Euclidean residuals | no |

RKCFC is therefore robust to non-Gaussian principal component scores and to a
covariance operator without finite fourth moments. It is not distribution-free in
general, and it is not protected against gross outliers in the observed values.

## Convergence

The iteration targets `J(c) = sum_i sum_j (Y_ij - Yhat_ij,(c_i))^2`, the total
squared reconstruction error at the current assignment. The reassignment step
minimises `J` exactly and is monotone by construction. The re-estimation step is
not, because the Kendall estimator does not minimise a sum of squared
reconstruction errors, and the smoothing bias does not vanish in `J`. Monotone
descent of the whole iteration is not guaranteed and neither is convergence of
the partition, and neither was established for the original KCFC.

In practice, run several starts and keep the partition with the smallest `J`.
Read `conv` together with `iter`, because `iter == maxIter` means the cap bound
the fit while `iter < maxIter` with `conv = FALSE` means a class fell below ten
curves. Those need opposite remedies.

## Also included

The metric, preprocessing and simulation machinery used by the accompanying
manuscript ships here, so that the version producing the published numbers is
the version that installs. The scripts of the simulation study and of the
real-data analysis are in the companion repository
[RKCFC_Simulation](https://github.com/BingbingFan/RKCFC_Simulation).

- `rre()`, the relative reconstruction error, where 0 is perfect and 1 is no
  better than the mean of the observations. It is unit-free, so simulation and
  application are comparable, which a raw mean squared error is not.
- `robust_normalize()` and `back_transform()`, median/IQR rescaling and its
  inverse for plotting in original units.
- `gen_scores()`, six score distributions ordered by how many moments exist,
  each calibrated to a common interquartile range so that tail behaviour is not
  confounded with signal-to-noise.
- `select_K()`, `ch_index()`, and `silhouette_index()`, score-space
  cluster-validity criteria, with a robust Mahalanobis whitening.
- `spatial_median()`, `sscm()`, `robust_scatter()`, `principal_angles()`.

## Version note

Versions before 0.2.0 contained two defects that changed reported results.

- `KCFC_init()` called the non-robust PACE estimator unconditionally, so
  `robust = TRUE` was started from a non-robust estimator. Since the iteration
  converges to a local optimum, the start selects the basin.
- `silhouette_score()` omitted the covariance bandwidth in its robust branch, so
  the class mean functions behind the coefficient were smoothed differently from
  the model being scored. `robust = FALSE` was unaffected.

Do not compare results across the boundary, and do not reproduce published
numbers with 0.1.x.

## Related

The multivariate extension, for clustering several longitudinal markers jointly,
is distributed separately as `MRKCFC`.

## Attribution

The Kendall FPCA runtime is derived from KFPCA version 2.0 by Rou Zhong and
Jingxiao Zhang (<https://github.com/BingbingFan/KFPCA>), incorporated under
GPL (>= 3) and identified in the source. RKCFC-specific adaptations are by
Bingbing Fan.

## References

Chiou, J.-M. and Li, P.-L. (2007). Functional clustering and identifying
substructures of longitudinal data. *JRSS-B* **69**(4), 679–699.

Zhong, R., Liu, S., Li, H. and Zhang, J. (2022). Robust functional principal
component analysis for non-Gaussian longitudinal data. *Journal of Multivariate
Analysis* **189**, 104864.

Yao, F., Müller, H.-G. and Wang, J.-L. (2005). Functional data analysis for
sparse longitudinal data. *JASA* **100**(470), 577–590.

## License

GPL (>= 3).
