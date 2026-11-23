# coresynth

`coresynth` は合成コントロール法（SCM）とその派生手法（SDID / GSC / MC / TASC / SI）を `RcppArmadillo` で実装した R パッケージ。全手法を統一の Formula インターフェース（`scm_fit()`）で提供し、最適化（QP・SVD・カルマンフィルタ）は C++ 実装。

## ビルド・テスト・検証

- **R 実行環境**: `C:\Program Files\R\R-4.6.2\bin\R.exe` / `Rscript.exe`（`bin\x64` ではなく `bin\` 直下）
- **ビルド・インストール・テスト**: `build-install` スキルを実行（Windows 特有の DLL ロック解除・`R CMD INSTALL` のサイレントロールバック回避手順込み）。テストは `tests/testthat/test-integration.R`（オブジェクトモデルのみ `test-object-model.R`）に集約
- **速度回帰チェック**: 計算経路や最適化に触れる更新の後は、関連する既存手法・条件（単一フィット、推論/プラセボ、予測子有無、パネル規模など）を多角的に実測し、精度とあわせて速度低下がないことを確認する（計測規則・基準は `.agents/docs/performance.md`）
- **pkgdown 検証**: 新しい export 関数・S3 method 追加時はコミット前に `pkgdown-check` スキルで `_pkgdown.yml` を検証する（漏れがあるとサイトデプロイ全体が停止）
- **CRAN チェック**: `Rscript.exe -e "devtools::check(args = c('--no-manual', '--run-donttest'))"`

## ドキュメント構成（共通仕様）

- **パス連動ルール**: `.agents/rules/*.md`（Antigravity 用）と `.claude/rules/*.md`（Claude Code 用）に同一内容を同期配置。該当ファイルを編集・レビューする前に必ず対応ルールに従う
- **共通リファレンス**: `.agents/docs/`
  - `theory-crosscheck.md`: 論文照合一覧・他パッケージ比較の基準値・被覆率検証値
  - `performance.md`: ベンチマーク基準・速度回帰の検証手順
  - `backlog.md`: 未実装拡張・破壊的変更候補
- **論文原典**: `papers/`（gitignore）

## 全手法横断仕様

- **行列次元規約**: 全手法 `T × N`（行 = 時間、列 = ユニット）
- **Staggered Adoption 共通加重（Clarke et al. 2024）**:
  $$\hat\tau = \frac{\sum_g N_{tr,g} \cdot T_{post,g} \cdot \hat\tau_g}{\sum_g N_{tr,g} \cdot T_{post,g}}$$
  - **コントロール群**: `"never_treated"`（既定）= 統制ユニットのみ。`"clean"`（非推奨・警告）= never-treated ＋ 将来採用者（事後窓をパネル末尾まで平均するため将来採用者が反実仮想を汚染する）。再フィット型の不足は `.check_control_group()` がエラー（SCM/GSC/SI は 2 未満、SDID は 0）。**MC/TASC は補完型**のため単一フィットの観測マスクで表現し `control_group` を持たないが、集計は同じ Clarke 加重を使う
- **推論関数の対応表**:

| 手法 | sharp 推論 | staggered 推論 |
|---|---|---|
| SCM | `mspe_ratio_pval()` / `scm_rank_pval()` | `scm_inference()`（jackknife 既定 / wild bootstrap） |
| SDID | `sdid_inference()`（4 method）/ `sdid_event_study()` | `sdid_inference()` のみ |
| GSC | `gsc_boot()` / `gsc_inference()` | `gsc_inference()` |
| MC | — | — |
| TASC | —（本パッケージで唯一、推論ツールなし） | — |
| SI | `si_inference()`（multi-arm 含む） | 同左 |

- `conformal_inference()` は sharp fit 限定（staggered / multi-arm / TASC はエラー）。`jackknife_global` は staggered 専用
- `sdid_inference()`・`gsc_inference()`・`si_inference()` は共通ドライバ `.unit_inference()` で**処置ユニットも再抽出・削除**する（`N_tr >= 2` 必須、staggered の jackknife は全コホート `N_tr,g >= 2` 必須）

## 環境・実装の絶対制約（落とし穴）

- **`arma::solve_opts::fast` 使用禁止**: Windows + Rtools45 でクラッシュする
- **OpenMP 並列化**: メインスレッドで `cs_num_threads()`（`src/threads.h`）を呼び `num_threads(...)` で上限を付ける（CRAN check 中は最大 2 スレッド、`options(coresynth.threads=)` 対応）。並列領域内は必ず `try/catch` で囲み（例外脱出で `std::terminate`）、R の状態アクセスや Rcpp 型（`Nullable` 等）を使わない
- **R / C++ 境界**: Armadillo `vec` は R 側で N×1 行列になる → R 層で `as.numeric()` 正規化
- **文字列結合の zero-length 罠**: `paste0("Donor ", integer(0))` は長さ 1 の `"Donor "` を返す → `sprintf()` を使う
- **非 ASCII 文字プロット禁止**: plot のタイトル・サブタイトル等にギリシャ文字や Unicode マイナスを入れない（非 UTF-8 ロケールの `R CMD check` でクラッシュ）。ASCII ハイフンと綴り（`lambda`/`omega`）を使う
- **リスト参照は `[[ ]]`**: `cohort_fits[[j]]` や fit 直下の要素に `$` を使うと部分マッチで別フィールドを吸う（例: `unit` → `unit_weights`、`F` → `force`、`se` → `separate_args`）
- **Formula の項制約**: `scm_fit()` は `all.vars()` で列名を取るため、変換式（`log(y)` 等）を黙って捨てないよう入口の `is.name()` 検査で拒否する

## コーディング・Git・リリース規約

- **命名規約**: C++ 関数 `snake_case_cpp`、R 内部ラッパー `fit_<method>_cpp()`、公開 API `scm_fit()`。ggplot `aes()` 変数は `utils::globalVariables()` に登録
- **コメント規約**: 論文節番号は `S.`（`§` 不使用）。開発履歴メモ（いつ・なぜ変えたか）はコードにもルールにも書かない（履歴はコミットと `NEWS.md` へ）
- **DESCRIPTION**: `Config/roxygen2/version` を維持（`devtools::document()` で `RoxygenNote:` が消えるのは roxygen2 8.0.0+ の仕様なので復元しない）
- **コミット・プッシュ**: 明示的な指示があるまで実行しない（作業後は変更を残して報告のみ）。Conventional Commits（英語・命令形・本文に why・AI 署名/Co-Authored-By なし）。`gh` CLI は `yo5uke` で認証済み。CRAN 提出用 `.tar.gz` 作成・アップロードはユーザー担当

## スコープ外（実装しない方針）

- multivariate outcome の同時推定、SI の time-varying covariates、Callaway & Sant'Anna 型 group-time CATT、Hirshberg & Klosin 型 conformal CI
