---
paths:
  - "R/gsc.R"
  - "src/gsc.cpp"
---

# GSC（Xu 2017）

$$Y_{it} = \delta_{it} D_{it} + x_{it}'\beta + \lambda_i' f_t + \varepsilon_{it}$$

共変量なしは truncated SVD 3 ステップ、共変量ありは EM ループ（E-step: SVD、M-step: リッジ OLS）。推論は `gsc_boot()`（sharp のみ、H0 下の parametric bootstrap）と `gsc_inference(method="bootstrap"/"jackknife"/"jackknife_global")`（sharp + staggered、非パラメトリック）。

**Staggered 設計**: `Y_co_g` は全 T 期間を使う（future-treated の処置後汚染を SDID と同様に許容する設計）。共変量は全ユニット T×N×p 配列を per-cohort にスライスして EM に渡す。

**因子正規化の注意**: 実装は `F = U_k · diag(D_k)` で Xu (2017) の識別制約（`F'F = T·I_r`）を満たさないが、反事実推定は F の列空間への射影であり正規化定数に依存しないため問題ない。
