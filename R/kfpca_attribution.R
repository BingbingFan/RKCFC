# KFPCA-derived internal runtime
#
# The following internal source files are derived from KFPCA version 2.0:
# basicfun.R, kernfun.R, GCVbw1D.R, GCVbw2D.R, GetGCVbw1D.R,
# GetGCVbw2D.R, GetMatList.R, MeanEst.R, MeanObsGrid.R, FPCscoreLSE.R,
# KFPCA.R, and predict.KFPCA.R.
#
# Upstream KFPCA authors: Rou Zhong and Jingxiao Zhang.
# Upstream source: https://github.com/BingbingFan/KFPCA
# Upstream license: GPL (>= 3).
#
# RKCFC integration and modifications by Bingbing Fan:
# - incorporated the runtime closure into the RKCFC namespace;
# - replaced direct generalized-inverse score calculations in FPCscoreLSE()
#   and predict.KFPCA() with RKCFC's solve()-first covariance solver;
# - retained these functions as internal implementation details of RKCFC.

# NOTE ON THE UPSTREAM EXAMPLES. The roxygen blocks in the KFPCA-derived files
# retain their upstream @examples, which call GenDataKL(). That generator is NOT
# vendored here: only the estimation runtime is. Those blocks are marked @noRd
# and so are never rendered or executed; they are kept because they document how
# the upstream authors intended each function to be called. KFPCA() and
# predict.KFPCA() are the two exported entry points and carry examples that run
# against pbc2 instead.
