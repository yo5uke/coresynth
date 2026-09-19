---
paths:
  - "R/mc.R"
  - "src/mc.cpp"
---

# MC（Athey et al. 2021）

$$\min_L \frac{1}{2} \|\mathcal{O} \odot (Y - L)\|_F^2 + \lambda_L \|L\|_*$$

反復 SVD ソフト閾値化（Soft-Impute）。λ_L は `|O|` で正規化されない（デフォルト = 0.01 × σ_max(Y)）。推論関数は無く、`conformal_inference()` のみ対応。

**観測マスク**: `O` は処置後セルと NA セルの**両方**を 0 にする（`O[is.na(Y)] <- 0`）。NA を観測値 0 として soft-impute に入れると ATT が大きく歪む。

**staggered はマスク経由**（TASC と共通）: `O` が各処置ユニットを自身の採用時点から 0 にするため、単一の分解がそのまま staggered fit になる。コホート別の再フィットはしない — `.impute_cohort_fits()`（`utils.R`）が単一 fit の `gap` を採用時点で分割してコホート表を作るだけ。したがって全コホートが同じドナー集合（never-treated ＋ 他処置ユニットの採用前セル）を共有し、`control_group` 引数は持たない（渡すとエラー）。集計は他手法と同じ Clarke 加重。**ユニット等加重にしてはいけない** — コホートの事後長が異なると同一データで他手法と別の推定量になる。

**ATT 過大推定の注意**: 小規模パネル（N_co ≈ 20 程度）では核ノルム正則化が強く効き、ATT が真値を上回りやすい。Mazumder et al. (2010) で理論的に知られた性質で、実装の不具合ではない。
