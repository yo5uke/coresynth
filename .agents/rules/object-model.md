---
trigger: glob
globs: "R/scm_fit.R,R/accessors.R,R/broom.R,R/unit_inference.R,R/export.R"
---

# fit オブジェクト・アクセサ・共通推論ドライバの契約

## クラス構造と S3 dispatch

- クラス付与は `new_coresynth()`（`R/scm_fit.R`）に一元化する。multi-arm SI では `"coresynth_multiarm"` を `"coresynth_staggered"` や `"coresynth_si"` より前に置き、`NextMethod()` で段階委譲できるようにする
- `staggered` / `multi_arm` の論理フィールドは後方互換のため残すが、ディスパッチは S3 クラスで行う
- `sdid_inference` は `c("sdid_inference", "coresynth_inference")` の二重クラスとし、`tidy()` / `glance()` は単一実装で全推論クラスを扱う

## staggered 集計と返却値の不変条件

- **集計の一元化**: コホート集計は全手法 `.aggregate_cohorts()`（`R/utils.R`）を通す。各要素は `cohort`・`n_treated`・`T_pre`・`T_post`・`estimate`・`weight`（`= N_tr,g * T_post,g`）を持ち、`weight` は表側で総和 1 に正規化する。SCM が `mass` / `unit_masses` を持つ場合は別名を作らず `.aggregate_cohorts()` 経由で `mass` を書き戻す
- **ブロックラベルの一元化**: `.block_labels()`（`R/utils.R`）を使う。`pool_cohorts = FALSE` では同一採用時点のブロックが複数並び `tidy()$term` や `augment()$.unit` が衝突するため、その場合だけ `unit_<id>` に切り替える。`cohort_estimates` の `unit` 列は、どのブロックもユニット名を持たない手法では列ごと落とす
- **MC / TASC の staggered 例外**: 単一フィットのため `Y_hat`・`gap` を保持するが、単一系列に潰すと採用時点の異なるコホートが混ざるため、`synthetic_outcomes.coresynth_staggered()` が `NULL` を返して `plot()` を明示エラーにする（系列確認は `augment()` を使う）
- **staggered `print` / `glance` の `T_pre`**: fit 直下の `T_pre` は最早採用コホートの値に過ぎないため staggered `print()` では出さない（コホート別は `tidy()`）
- **SDID の `.resid` 事後平均 ≠ `estimate`**: SDID の ATT は λ 加重の事前期間を引いた DiD であり、描画される合成系列とのギャップ平均とは λ が一様からずれる分だけ一致しない。これは推定量固有の性質（バグではない）

## アクセサとエクスポート

- 手法別のドナー行列フィールド名（SCM/SDID: `Y_co_pre`/`Y_co_post`、SI: `Y_pre_co`/`Y_post_co`、GSC/MC: `Y_co_all`）は `accessors.R` だけが知り、`donor_outcomes.coresynth()` は保存済み旧オブジェクト向けに全レイアウト探索のフォールバックを保つ
- `cohort_fits[[j]]` や fit 直下からフィールドを読むときは**必ず `[[ ]]`** を使う（`$` は `unit` → `unit_weights`、`F` → `force`、`se` → `separate_args` のように部分マッチする）

## 共通推論ドライバ `.unit_inference()`（`R/unit_inference.R`）

- GSC・SI・SDID が共有する。各手法は `estimate_fn(ctx, co, tr)`（重みを固定する SDID jackknife は別に `jack_fn`）を渡し、fit と同値再現（`1e-10` 許容）を保証する
- **処置ユニットも必ず再抽出・削除する**（対照だけだと処置側の事前ノイズ・事後ショックが SE から抜ける）。`N_tr = 1`（multi-arm はいずれかの arm が 1）はエラー、staggered の jackknife は全コホート `N_tr,g >= 2` 必須（Clarke et al. 2024 Algorithm 3）
- **bootstrap の抽出規則の違い**: GSC・SI は処置（arm 層別）と対照を層別抽出、SDID は `ctx$pooled` で全ユニットをまとめて抽出し処置 0 または対照 0 なら引き直す（Arkhangelsky Algorithm 2）。いずれも引いた対照を全コホートで共有し、引いた処置からコホートと Clarke 重みを組み直す
- **jackknife**: ユニットごとに 1 回だけ再推定し、per-cohort（delta-method）と `jackknife_global` の両方へ使い回す。staggered では fit の `control_group` を引き継ぐ（記録のない旧 fit は `.fit_control_group()` で判定）
