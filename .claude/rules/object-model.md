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
- **MC/TASC は例外**で `Y_hat`・`gap` を保持する（単一フィットなので per-unit 反実仮想が実在する）。ただし単一系列に潰すと採用時点の違うコホートを平均してしまうため、`synthetic_outcomes.coresynth_staggered` が NULL を返して `plot()` を明示エラーに倒す。系列は `augment()` のコホート別 long df を使う
- 集計は全手法 `.aggregate_cohorts()`（`utils.R`）に一元化。各 cohort_fit が `cohort`/`n_treated`/`T_pre`/`T_post`/`estimate`/`weight`（= N_tr,g × T_post,g）を持つことが前提で、`weight` は表側で総和 1 に正規化される。SCM は加えて `mass`（ブロック質量）/`unit_masses` を持つ — `unit_weights` 経由で呼ばれたときは `.aggregate_cohorts()` も `mass` を書き戻すので**同じ量の別名を作らない**
- **ブロックのラベルは `.block_labels()`（`utils.R`）に一元化**。`pool_cohorts = FALSE` では同一採用時点のブロックが複数あるため、コホート名だと `tidy()$term` が重複し `augment()$.unit` が**別ユニットの系列を 1 本に潰す**。その場合だけ `unit_<id>` に切り替える。`cohort_estimates` には `unit` 列が付く（1 ユニットのブロックのみ非 NA、どのブロックも名乗らなければ列ごと落とす → 他手法の表はバイト不変）。`cohort_fits[[j]]` の `unit` は `[[ ]]` で読む — SDID が `unit_weights` を持つので `$` が部分マッチする
- `print`/`summary` は `coresynth_staggered` メソッドが `(staggered)` ラベルとコホート表を出す。fit 直下の `T_pre` は**最早採用コホートの値**でしかないため staggered の print では出さない（`glance()` の `T_pre`/`T_post` は同じ理由でパネル水準の値であり、コホート別は `tidy()` を使う）

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
- **SDID だけ `.resid` の事後平均が `estimate` と一致しない**（sharp・staggered とも）。SDID の ATT は λ 加重のプリ期間を引いた DiD なので、描画される synthetic 系列とのギャップ平均とは λ が平坦からずれる分だけ違う。他の 5 手法は反実仮想がすべてなので一致する。実測でも 1.874 vs 1.808（sharp）・1.649 vs 1.783（staggered の遅いコホート）。**推定量固有の性質でありバグではない**ので、`estimators.Rmd` の SDID 節に明記するにとどめている

## staggered 推論の共通ロジック（per-cohort 集約）

- **staggered jackknife（per-cohort）**: per-cohort LOO + delta-method 分散合成（`sum((w_g/W)^2 * var_g)`）。コホート間相関を無視
- **jackknife_global**: 対照も処置も含むユニーク集合を横断 LOO。同一ユニットが複数コホートに属す場合（対照、`"clean"` の将来採用者）は全コホートから同時除外 → 全コホート再推定 → 加重平均（落とした処置のコホートは重みも減る）。コホート間相関を捕捉する。`V̂ = (N-1)/N × Σ_i (ATT^(-i) - mean)^2`（N = 落としたユニット数）。staggered 専用（sharp・sharp multi-arm では不可）
- **GSC・SI・SDID は共通ドライバ `.unit_inference()`（`R/unit_inference.R`）**（1.1.0.9000）: 手法側は「ある対照集合・処置集合（重複可）で推定量を当て直す」関数 `estimate_fn(ctx, co, tr)` を持つ（`.gsc_unit_estimate()`/`.si_unit_estimate()`/`.sdid_unit_estimate()`、fit の再現は `expect_equal(..., tolerance = 1e-10)` でテスト済み）。jackknife が重みを固定する SDID は別に `jack_fn`（`.sdid_jack_estimate()`）を渡す。**処置も再抽出・削除の対象**で N_tr=1 は明示エラー。**staggered の jackknife はどのコホートにも処置 2 以上が必要**（Clarke et al. 2024 Algorithm 3、3 手法共通でエラー）。bootstrap は GSC・SI が処置を層（multi-arm のアーム）内で引き対照を別に引く層別抽出（fect と同じ）、SDID は `ctx$pooled` で全ユニットをまとめて引いて処置 0・対照 2 未満を引き直す（Arkhangelsky Algorithm 2 の字面）。どちらもコホートはプールし、引いた対照を全コホートで共有、引いた処置からコホートと Clarke 重みを組み直す。jackknife は対照も処置も落とし、ユニットごとの再推定 1 回を global と per-cohort の両方に使い回す。fit は `control_group` を記録し、記録のない旧 fit は `.fit_control_group()` が cohort_fits から判定する
