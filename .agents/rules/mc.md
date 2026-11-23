---
trigger: glob
globs: "R/mc.R,src/mc.cpp"
---

# MC（Athey et al. 2021, JASA; Mazumder et al. 2010, JMLR）

$$\min_{L, a, b} \frac{1}{2} \|\mathcal{O} \odot (Y - L - a \mathbf{1}_N^\top - \mathbf{1}_T b^\top)\|_F^2 + \lambda_L \|L\|_*$$

- **パラメータ化と観測マスク**: `lambda_L` は `|O|` で正規化しない（Athey et al. 2021 Eq. 4.3 の `lambda * |O| / 2` を直接表す）。観測マスク `O` は処置後セルと `NA` セルの**両方**を `0` にする（`.conformal_refit.coresynth_mc()` も同様）
- **staggered と識別不能セル**: staggered は `O` が各処置を採用時点から `0` にすることで単一分解がそのまま対応し、`.impute_cohort_fits()` が分割して Clarke 加重で集計する（`control_group` 引数は持たない。ユニット等加重禁止）。未処置観測ゼロの期間および未処置期間ゼロのユニットは推定値が 0 に張り付いて水準バイアスを生むため、`.drop_unidentified()`（`R/utils.R`）で警告つき除外し、切り詰め後のパネルから組み直す
- **二元固定効果（`fixedeff = TRUE` 既定、`soft_impute_fe_cpp`）**: Athey et al. 2021 Eq. 4.3（および MCPanel 既定）に従い非罰則 FE を推定する（FE を外すとアウトカム水準が `L` に吸われ核ノルム縮退で ATT が過大になる）。マスク下で直交しないため 1 反復あたり 5 スイープのバックフィットと SVD を同一ループ内で交互実行する。Choi & Yuan (2024) の脱バイアスは MNAR 用であり FE 単独に対し改善しないため不採用
- **収束判定は `L` 単独の相対変化**: 完成行列 `M = L + FE` で判定すると水準を含む分だけ早く停止し精度が落ちるため、`L` のみで判定する
- **λ の CV（`.mc_cv_lambda()`）**:
  - 罰則は `L` のみに掛かるため、パス上端 `.mc_lambda_max()` は `fixedeff = TRUE` なら「`L = 0` の TWFE 残差」の最大特異値、FALSE なら生パネルの最大特異値を使う（生パネルを使うと高水準データでランク 0 に退化し水準不変性が壊れる）。フォールバック `0.01 * sigma_max` も同じ上端を使う
  - 観測率を揃えるため `|O_k| = floor(|O|^2 / (N*T))` で fold を引き、**自前 RNG ストリームで `.Random.seed` を保存・復元**する。warm-start は `soft_impute_cpp` のゼロ初期化定義とずれるため不採用。fold が成立しない小パネルではフォールバックし `fit$lambda_cv = FALSE` を記録
- **推論ツール**: `conformal_inference()` のみ対応
