# coresynth: Fast and Unified Synthetic Control Methods

A unified 'Formula' interface to the Synthetic Control Method (SCM) and
related panel-data causal inference estimators: Synthetic
Difference-in-Differences (SDID), Generalized Synthetic Control (GSC),
Matrix Completion (MC), Time-Aware Synthetic Control (TASC), and
Synthetic Interventions (SI), together with an experimental-design
variant. Computational bottlenecks (quadratic programming, singular
value decomposition, and Kalman filtering) are implemented in 'C++' via
'RcppArmadillo'. Methods are described in Abadie, Diamond and
Hainmueller (2010)
[doi:10.1198/jasa.2009.ap08746](https://doi.org/10.1198/jasa.2009.ap08746)
, Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021)
[doi:10.1257/aer.20190159](https://doi.org/10.1257/aer.20190159) , Xu
(2017) [doi:10.1017/pan.2016.2](https://doi.org/10.1017/pan.2016.2) ,
Athey, Bayati, Doudchenko, Imbens and Khosravi (2021)
[doi:10.1080/01621459.2021.1891924](https://doi.org/10.1080/01621459.2021.1891924)
, and Agarwal, Shah and Shen (2025)
[doi:10.1287/opre.2025.1590](https://doi.org/10.1287/opre.2025.1590) .

## See also

Useful links:

- <https://github.com/yo5uke/coresynth>

- <https://yo5uke.com/coresynth/>

- Report bugs at <https://github.com/yo5uke/coresynth/issues>

## Author

**Maintainer**: Yosuke Abe <yosuke.abe0507@gmail.com>

Authors:

- Yosuke Abe <yosuke.abe0507@gmail.com>
