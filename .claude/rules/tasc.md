---
paths:
  - "R/tasc.R"
  - "src/tasc.cpp"
---

# TASC（Rho et al. 2026）

$$z_{t+1} = A z_t + C + \eta_t, \quad y_t = W z_t + \varepsilon_t$$

EM アルゴリズム（E-step: カルマンフィルタ + RTS スムーザ、M-step: OLS）。Joseph form で数値安定化。**推論関数・`conformal_inference()` とも未対応**（本パッケージで推論ツールを持たない唯一の手法）。

**M-step（W・R）の注意**: 処置ユニットは事後期間が NA のため、OLS 更新は **per-unit ループ**で観測済み時点のみを使う必要がある（Shumway-Stoffer missing-data EM 準拠）。全 T 期間を分母に使うと処置ユニットの W が T_pre/T 倍にバイアスされ ATT を過大推定する。

**`Y_hat` の規約**: TASC だけ `Y_hat` が T×N の**全ユニット**行列（GSC は `Y_tr_hat`、MC は処置列のみを返す）。**消費側が処置列を選ぶ**のが正しい規約 — `synthetic_outcomes()`・`export_json()` は `idx_tr` 列だけを使う。全ユニット平均を使うと plot/augment/export の表示が壊れる（ATT 自体は per-unit gap 由来のため影響しない）。

**staggered フラグを持たない**: `fit_tasc_cpp` は Y_hat を常に全 T×N で返す設計のため、多時点採用でも cohort 集計をせず `staggered = FALSE` のまま動作する。
