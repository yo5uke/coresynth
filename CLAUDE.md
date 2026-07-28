# coresynth — CLAUDE.md

本ファイルは Claude Code が本プロジェクトで作業する際のコンテキストガイド。詳細な理論的根拠・アルゴリズム決定事項は `.claude/rules/*.md`（該当ファイル編集時に自動読込）と `.claude/docs/*.md`（随時参照）に分離。

---

## プロジェクト概要

**coresynth** は合成コントロール法（SCM）とその派生手法を RcppArmadillo で高速実装した R パッケージ。

- すべての推定手法に **統一的な Formula インターフェース** を提供
- ボトルネックとなる最適化（QP、SVD、カルマンフィルタ）を全て **C++ で実装**
- Synth パッケージ比の実測速度は `.claude/docs/performance.md` 参照（開発版ソルバ: N_co=40〜100 で単一フィット ~2,100〜2,700 倍・in-space プラセボ ~25,000〜26,000 倍）

### 対象手法

論文・数式・アルゴリズム決定事項は各手法の詳細ファイル参照（該当 R/C++ ファイルを開くと自動読込）。

| 手法 | 詳細 |
|---|---|
| SCM | `.claude/rules/scm.md` |
| SDID | `.claude/rules/sdid.md` |
| GSC | `.claude/rules/gsc.md` |
| MC | `.claude/rules/mc.md` |
| TASC | `.claude/rules/tasc.md` |
| SI | `.claude/rules/si.md` |
| SCM-Design | `.claude/rules/scm-design.md` |

横断的推論として **Conformal Inference**（`R/conformal.R`）が `scm`/`sdid`/`gsc`/`mc`/`si` の sharp fit に対応（詳細は `.claude/rules/conformal.md`）。

横断的ルール: `.claude/rules/object-model.md`（fit オブジェクトのクラス構造・アクセサ・staggered 推論の共通ロジック）、`.claude/rules/plot.md`（plot 内部規約）、`.claude/rules/input-validation.md`（入力検証）。随時参照: `.claude/docs/theory-crosscheck.md`（全手法横断の論文照合表）、`.claude/docs/performance.md`（ベンチマーク）。

論文テキスト・数式整理 MD は `papers/`（.gitignore 管理）に格納。

### Staggered Adoption 共通仕様（Clarke et al. 2023）

$$\hat\tau = \frac{\sum_g N_{tr,g} \cdot T_{post,g} \cdot \hat\tau_g}{\sum_g N_{tr,g} \cdot T_{post,g}}$$

- **コホート g**: 採用タイミングが同じ処置ユニット群
- **コントロール群**: `"clean"`（デフォルト）= never-treated ＋ 将来採用者（$T_{adopt}[j] > g$）、`"never_treated"` = 統制ユニットのみ
- **返却構造**: `staggered = TRUE`、`cohort_estimates`（data.frame）、`cohort_fits`（list）
- **cohort_fits[[g]]** の共通フィールド: `estimate`, `weight`, `idx_co`, `idx_tr`, `T_pre`, `T_post`（inference resampling に利用）。SCM は加えて `Y_treat_mat`（T × n_treated_g、per-unit 処置系列 — wild bootstrap 用）と `alpha`（fixedeff の intercept shift、なければ 0）

---

## ファイル構成

```
coresynth/
├── R/
│   ├── scm_fit.R          Formula インターフェース・dispatch・new_coresynth()
│   ├── scm.R              SCM・mspe_ratio_pval()・placebo_in_time()・loo_donors()
│   ├── sdid.R             SDID・sdid_inference()
│   ├── gsc.R              GSC・gsc_boot()・gsc_inference()
│   ├── mc.R               MC ラッパー
│   ├── tasc.R             TASC ラッパー
│   ├── si.R               SI・si_inference()・multi-arm 実装
│   ├── scm_design.R       SCM-Design（Abadie & Zhao 2026）
│   ├── conformal.R        conformal_inference()（CWZ 2021）
│   ├── utils.R            panel_to_matrices()・panel_to_tensor()・pred()・
│   │                      build_predictor_matrices()・.check_panel_complete()
│   ├── accessors.R        treated_outcomes()・synthetic_outcomes()・donor_outcomes()
│   ├── broom.R            tidy() / glance() / augment()
│   ├── plot.R             plot.coresynth()・plot.scm_placebo()（ggplot2）・plot_data()
│   ├── export.R           export_json()
│   └── coresynth-package.R useDynLib 宣言
├── src/
│   ├── optim.cpp          proj_simplex, solve_simplex_qp（FISTA）
│   ├── scm.cpp            scm_inner_weights_cpp, scm_weights_cpp
│   ├── sdid.cpp           sdid_unit/time_weights_cpp, sdid_estimate_cpp
│   ├── gsc.cpp            gsc_ife_cpp（EM ループ・共変量対応）
│   ├── mc.cpp             soft_impute_cpp
│   ├── tasc.cpp           kalman_smoother_cpp
│   ├── si.cpp             tensor_unfold_cpp, si_pcr_cpp
│   └── inference.cpp      sdid_placebo_cpp, scm_placebo_cpp（OpenMP 並列化済み）
├── tests/testthat/
│   ├── test-integration.R   統合テスト（大半はここに集約）
│   └── test-object-model.R  クラスタグ・アクセサ・S3 dispatch
├── README.Rmd             ソース（knitr 生成物が README.md — README.md は直接編集しない）
├── .claude/rules/, .claude/docs/, .claude/skills/  Claude Code 向け詳細ドキュメント・スキル
└── sandbox/, papers/, plan/  gitignore（ベンチマーク・検証スクリプト / 論文 MD / 計画メモ）
```

---

## 推論関数の対応表

引数・返却値の詳細は各関数の roxygen ドキュメント（`?関数名`）参照。ここには複数手法・複数関数を横断しないと分からない対応関係だけをまとめる。

| 手法 | sharp 推論 | staggered 推論 | `conformal_inference()` |
|---|---|---|---|
| SCM | `mspe_ratio_pval()`（順列検定） | `scm_inference()`（wild bootstrap） | ✅ |
| SDID | `sdid_inference()`（4 method 全て） | 同左 | ✅ |
| GSC | `gsc_boot()` / `gsc_inference()` | `gsc_inference()` | ✅ |
| MC | — | — | ✅ |
| TASC | — | — | ❌（推論ツール自体なし） |
| SI | `si_inference()`（multi-arm 含む） | 同左 | ✅ |

`jackknife_global` は staggered 専用（sharp・sharp multi-arm では不可）。`conformal_inference()` は sharp fit 限定（staggered/multi-arm/tasc はエラーで上表の手法別関数へ誘導）。fit オブジェクトのクラス構造・アクセサ（`treated_outcomes()`/`synthetic_outcomes()`/`donor_outcomes()` 等）は `.claude/rules/object-model.md`、`plot()`/`plot_data()` の内部規約は `.claude/rules/plot.md` 参照。

---

## 実装上の落とし穴（環境・R・C++）

- **`arma::solve_opts::fast` は使用禁止** — 本環境（Windows + Rtools45）でクラッシュ
- Armadillo `vec` は R 側で **N×1 行列**になる → R 層で `as.numeric()` 正規化
- `paste0("Donor ", integer(0))` は zero-length が `""` 扱いになり長さ 1 の `"Donor "` を返す → zero-length を保持する `sprintf()` を使う
- `R CMD INSTALL` 失敗時は旧版に**サイレントロールバック**される — 失敗に気づかず古い DLL でテストする罠。ビルド・テストは `build-install` skill を使う（Stop-Process 手順込みで実行される）
- 本マシンに gh CLI はない。PR 作成は `git credential fill` + GitHub REST API で行う
- **plot のタイトル・サブタイトルに非 ASCII 文字を入れない**（ギリシャ文字・Unicode マイナス記号 `−` U+2212 等）— UTF-8 でない ロケール（GitHub Actions macOS ランナー等）で grid のテキストメトリクス計算が `mbcsToSbcs` 変換に失敗し `R CMD check --run-donttest` がクラッシュする。マイナス記号は ASCII ハイフン `-`、ギリシャ文字は綴り（`lambda`/`omega`）で代替する。R ソース自体は ASCII のため `\uXXXX` エスケープで紛れ込みやすい

---

## ビルド手順

`build-install` skill が正確な手順（DLL ロック解除・コピー込み）を実行する。手動で行う場合:

```powershell
# R ファイルのみ変更後
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL .

# C++ 変更後（DLL ロック解除が必要）
Get-Process | Where-Object {$_.Name -match "Rterm"} | Stop-Process -Force
Start-Sleep -Seconds 1
Remove-Item src\*.o -Force -ErrorAction SilentlyContinue
& "C:\Program Files\R\R-4.6.1\bin\R.exe" CMD INSTALL .

# テスト実行（DLL コピーが必要）
Copy-Item "$env:USERPROFILE\AppData\Local\R\win-library\4.6\coresynth\libs\x64\coresynth.dll" `
  -Destination "src\coresynth.dll" -Force
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "devtools::test()"
```

> R は 4.6.1（2026-07 時点、`bin\x64` サブディレクトリではなく `bin\` 直下の実行ファイルを使用）。

---

## コミット

明示的な指示があるまでコミット・プッシュは行わない（作業後は変更を残して報告のみ）。

### 新規 export 追加時のチェックリスト

新しいエクスポート関数・S3 method（`plot.*`、`print.*` 等）を追加したら、`pkgdown-check` skill でコミット前に検証する（`_pkgdown.yml` の reference index にトピックが無いと GitHub Actions の pkgdown デプロイがサイト全体停止するため）。手動で行う場合の手順・検証コマンドは `.claude/skills/pkgdown-check/SKILL.md` 参照。

---

## コーディング規約

- C++ 関数名: `snake_case_cpp` サフィックス（例: `sdid_unit_weights_cpp`）
- R ラッパー: `fit_<method>_cpp()` が内部関数、`scm_fit()` が公開 API
- 行列の次元規約: `T × N`（行=時間、列=ユニット）統一
- テストは `tests/testthat/test-integration.R` に集約（オブジェクトモデルのみ `test-object-model.R`）
- ggplot の `aes()` 変数は `utils::globalVariables()` に登録（R CMD check NOTE 回避）
- コメント: 論文セクション参照は `§` でなく `S.`。開発履歴メモ（いつ・なぜ変更したか）はコメントや CLAUDE.md に書かず、コミットメッセージ・NEWS.md に書く

### roxygen2 規約（バージョン 8.0.0）

```
Config/roxygen2/markdown: TRUE   # DESCRIPTION に記載
RoxygenNote: 8.0.0
```

`devtools::document()` 後に `RoxygenNote:` が 7.x に書き戻された場合は即座に `8.0.0` に修正する。

---

## 既知の制限

推論関数の可否は上表参照。それ以外の制限:

- **SI の共変量**: 未対応（sharp/staggered とも。Agarwal et al. 2025 の原著が covariate なしの定式化）
- **SCM staggered**: `predictors = NULL` 限定（`covariates` partial-out のみ。エラー時に代替案を案内）
- **plot()**: staggered fit は `Y_synth = NULL` のため未対応（明示エラー）。cohort 別の系列は `augment()` の long df を使う

---

## スコープ外（意図的な制限）

応用論文で頻出するが本パッケージでは対応しない:
- multivariate outcome（同時推定）— Dube (2026) の "ウェイト流用" 用途では使えない
- SI の time-varying covariates — Agarwal et al. (2025) は covariate なしの定式化が原著
- Callaway & Sant'Anna (2021) 型 group-time CATT — 別フレームワーク（コホート ATT 集計は対応済み）
- Hirshberg & Klosin (2023) 型 conformal CI — 別枠

---

## その他

- コミット時、コミットメッセージにはClaudeのサインを記載しない。