---
paths:
  - "R/scm_fit.R"
  - "R/accessors.R"
  - "R/broom.R"
---

# fit オブジェクトのクラス構造

`new_coresynth()` コンストラクタ（`R/scm_fit.R`）が一元付与。sharp = `c("coresynth_<method>", "coresynth")`、staggered は先頭に `"coresynth_staggered"`、multi-arm SI はさらにその前に `"coresynth_multiarm"`（multi-arm メソッドが勝って `NextMethod()` で委譲できる順序）。`staggered`/`multi_arm` のリストフィールドは後方互換のため保持。print/summary/tidy/augment/`.conformal_refit` は S3 dispatch（`.conformal_refit` は estimator ごとに 1 メソッドの generic）。

**staggered fit の共通返却フィールド**:
- `staggered = TRUE`, `cohort_estimates` (data.frame), `cohort_fits` (list)
- `estimate`: 加重平均 ATT, `Y_treat`: T × N_tr, `Y_synth = NULL`, `gap = NULL`
- SDID のみ: `Y_mat` (T × N_all, residualized — inference の per-cohort resample 用)

**multi-arm fit（SI のみ、sharp）の返却フィールド**:
- `multi_arm = TRUE`, `arm_levels` (integer vector, 制御 arm 0 を除く), `arm_estimates` (named numeric)
- `arm_fits` (named list): per-arm に `weights` (N_co × N_tr_d), `Y_cf`, `Y_treat`, `Y_synth`, `gap`, `idx_tr`, `idx_co`, `T_pre`, `T_post`, `estimate`, `weight`
- `Y_pre_co`, `Y_post_co`, `idx_co` (inference 用), `staggered = FALSE`

**staggered × multi-arm fit（SI のみ）の返却フィールド**:
- `staggered = TRUE`, `multi_arm = TRUE`, `arm_levels`, `arm_estimates`（arm ごとのコホート加重平均）
- `cohort_arm_estimates`（data.frame: cohort・arm 両列）、`cohort_fits`（フラットリスト、各 (g,d) セルに cohort/arm フィールド）
- `Y_all`（T × N — inference 用）。`tidy()` は `term="cohort_g_arm_d"` 形式で返す

**`coresynth_inference` / `sdid_inference` 共通フィールド**:
- `estimate`, `se`, `p_value`, `ci_lower`, `ci_upper`
- `method`, `staggered` (logical), `n_controls`, `alternative`
- `boot_ests`（bootstrap 時: 数値ベクトル; jackknife 時: NULL）
- `sdid_inference` は `c("sdid_inference", "coresynth_inference")` 二重クラス。`tidy()`/`glance()` は単一実装で全推論クラスを扱う

## アクセサ規約

- **アクセサ**: 手法別のドナー行列フィールド名（SCM/SDID: `Y_co_pre`/`Y_co_post`、SI: `Y_pre_co`/`Y_post_co`、GSC/MC: `Y_co_all`）はアクセサだけが知る。旧バージョン保存オブジェクト向けに `donor_outcomes.coresynth` は全レイアウト探索のフォールバックを持つ
- **augment staggered**: SDID cohort_fits は Y_synth を持たないため `Y_all_parent %*% omega + omega0` で反実仮想を再構成（SCM/SI は Y_synth を直接保持）。multi-arm staggered は `.arm` 列を追加。MC/TASC は unit_weights を持たないため `glance()` の `n_controls` は `idx_co`/`idx_tr` fallback

## staggered 推論の共通ロジック（per-cohort 集約）

- **staggered jackknife（per-cohort）**: per-cohort LOO + delta-method 分散合成（`sum((w_g/W)^2 * var_g)`）。コホート間相関を無視
- **jackknife_global**: ユニーク control 集合を横断 LOO。同一 control が複数コホートに属す場合は全コホートから同時除外 → 全コホート再推定 → 加重平均。コホート間相関を正確に捕捉。`V̂ = (N-1)/N × Σ_i (ATT^(-i) - mean)^2`（N = unique control 数）。staggered 専用（sharp・sharp multi-arm では不可）
