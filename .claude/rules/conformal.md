---
paths:
  - "R/conformal.R"
---

# Conformal Inference（Chernozhukov, Wüthrich & Zhu 2021）

`scm`/`sdid`/`gsc`/`mc`/`si` の sharp fit に対応する横断的推論。**既知の制限**: sharp のみ（staggered/multi-arm/tasc は明示エラーで既存の `*_inference()` へ誘導）。

**設計上の要点**: 帰無 H0:τ=τ0 のもとで処置後を `Y1-τ0` に補完し、**全 T 期間で反実仮想を再推定**する（事前期間のみだと有限標本妥当性が崩れる、CWZ §2.2）。`.conformal_refit` は S3 generic で、手法ごとに「全 T 期間」に対応する再推定ロジックを持つ（SCM=全 T の simplex QP、SI=全 T の SVD、SDID=全 T のユニット重み＋concentrated intercept、GSC=EM を T_pre=T として実行、MC=soft-impute を全観測で実行）。共変量つき fit は残差化済み行列でそのまま動く。検定統計量は残差の moving-block（巡回シフト）置換で p 値を、test inversion グリッドで CI を作る。
