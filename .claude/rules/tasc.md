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

**staggered はマスク経由**（MC と共通）: `Y_obs` が各処置ユニットを自身の採用時点から NA にするため、単一の状態空間フィットがそのまま staggered fit になる。コホート別の再フィットはしない — `.impute_cohort_fits()`（`utils.R`）が単一 fit の `gap` を採用時点で分割してコホート表を作るだけ。したがって全コホートが同じドナー集合を共有し、`control_group` 引数は持たない（渡すとエラー）。集計は他手法と同じ Clarke 加重。**ユニット等加重にしてはいけない** — コホートの事後長が異なると同一データで他手法と別の推定量になる。SVD 初期化だけはグローバル `T_pre`（最早採用時点 - 1）を使うが、EM 本体が per-unit マスクで補正するため推定には効かない。
