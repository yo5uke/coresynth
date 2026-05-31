# coresynth: Fast and Unified Synthetic Control Methods

High-performance Synthetic Control Method (SCM) and related causal
inference estimators (SDID, GSC, MC, TASC, SI, plus an
experimental-design variant) with a unified Formula interface. All
computational bottlenecks (QP solving, SVD-based matrix completion,
Kalman filtering) are implemented in C++ via RcppArmadillo, providing up
to ~55x speed improvements over pure-R alternatives on typical problem
sizes.

## See also

Useful links:

- <https://github.com/yo5uke/coresynth>

- <https://yo5uke.com/coresynth/>

- Report bugs at <https://github.com/yo5uke/coresynth/issues>

## Author

**Maintainer**: Yosuke Abe <yosuke.abe0507@gmail.com>

Authors:

- Yosuke Abe <yosuke.abe0507@gmail.com>
