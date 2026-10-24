# coresynth

**coresynth** は合成コントロール法（SCM）とその派生手法（SDID/GSC/MC/TASC/SI）を RcppArmadillo で高速実装した R パッケージ。全手法を統一的な Formula インターフェース（`scm_fit()`）で提供し、最適化（QP・SVD・カルマンフィルタ）は C++ 実装。

手法ごとの論文根拠・アルゴリズム決定事項は `.claude/rules/<method>.md`（`scm`/`sdid`/`gsc`/`mc`/`tasc`/`si`/`scm-design`/`conformal`）にある。**該当手法の R/C++ を編集・レビューする前に必ず読む**。横断ルールは `object-model.md`（fit オブジェクトのクラス構造・アクセサ・staggered 推論の共通ロジック）、`plot.md`、`input-validation.md`、`versioning.md`。ここには**複数ファイル・複数手法を横断しないと分からないこと**だけを書く。

自動では読み込まれない参照資料: `.claude/docs/theory-crosscheck.md`（論文照合表・他パッケージとの基準値）、`.claude/docs/performance.md`（ベンチマーク）、`.claude/docs/backlog.md`（保留項目）。論文 PDF/MD は `papers/`（gitignore）。

---

## Staggered Adoption 共通仕様（Clarke et al. 2024）

$$\hat\tau = \frac{\sum_g N_{tr,g} \cdot T_{post,g} \cdot \hat\tau_g}{\sum_g N_{tr,g} \cdot T_{post,g}}$$

- **コントロール群**: `"never_treated"`（既定）= 統制ユニットのみ。`"clean"`（非推奨・警告、次メジャーで削除）= never-treated ＋ 将来採用者。**コホート効果は採用時点からパネル末尾まで平均するので、将来採用者は事後窓の途中で処置済みになり反実仮想を汚染する**（τ=3 で −0.24〜−0.59 のバイアスを実測）。Clarke Algorithm 1・BFR S.2.2 とも、末尾まで推定するなら never-treated に一致する。never-treated が 2 未満なら `.check_control_group()` が明示エラーで MC/TASC へ誘導する。**MC/TASC は `control_group` を持たない**
- **コホート分割する手法としない手法**: SCM/SDID/GSC/SI はコホートごとに再フィットする。**MC/TASC は補完型**なので再フィットせず、観測マスクが採用時点を表現し単一フィットがそのまま staggered 推定になる。**どちらも集計は同じ Clarke 加重**でなければならない — ユニット等加重にすると事後長の異なるコホートで他手法と別の推定量を返す

## 推論関数の対応表

| 手法 | sharp 推論 | staggered 推論 |
|---|---|---|
| SCM | `mspe_ratio_pval()` | `scm_inference()`（jackknife 既定 / wild bootstrap） |
| SDID | `sdid_inference()`（4 method 全て。bootstrap・jackknife は N_tr ≥ 2、1 なら placebo）/ `sdid_event_study()`（期間別効果＋クラスタ/ユニット bootstrap） | `sdid_inference()` のみ |
| GSC | `gsc_boot()`（Xu Algorithm 2、N_tr が少ないとき）/ `gsc_inference()`（N_tr ≥ 2） | `gsc_inference()` |
| MC | — | — |
| TASC | — | — （本パッケージで唯一、推論ツール自体なし） |
| SI | `si_inference()`（multi-arm 含む、N_tr ≥ 2。1 なら `conformal_inference()`） | 同左 |

`jackknife_global` は staggered 専用。`conformal_inference()` は sharp fit 限定（staggered/multi-arm/tasc はエラーで上表へ誘導）。`sdid_inference()`・`gsc_inference()`・`si_inference()` の bootstrap・jackknife は処置ユニットも再抽出・削除する（共通ドライバ `.unit_inference()`、`object-model.md`）。staggered の jackknife はどのコホートにも処置 2 以上が必要（Clarke Algorithm 3）。

---

## 実装上の落とし穴（環境・R・C++）

- **`arma::solve_opts::fast` は使用禁止** — 本環境（Windows + Rtools45）でクラッシュ
- Armadillo `vec` は R 側で **N×1 行列**になる → R 層で `as.numeric()` 正規化
- `paste0("Donor ", integer(0))` は zero-length が `""` 扱いになり長さ 1 の `"Donor "` を返す → zero-length を保持する `sprintf()` を使う
- `R CMD INSTALL` 失敗時は旧版に**サイレントロールバック**される。ビルド手順は下記
- **plot のタイトル・サブタイトルに非 ASCII 文字を入れない**（ギリシャ文字・Unicode マイナス記号等）— UTF-8 でないロケールで `R CMD check --run-donttest` がクラッシュする。マイナスは ASCII ハイフン、ギリシャ文字は綴り（`lambda`/`omega`）で代替する
- **`cohort_fits` から新しいフィールドを読むときは `[[ ]]`**（`$` は部分マッチする — SDID の `unit_weights` に `unit` が、GSC の `r_cv` に `r` が無言で吸われる）
- **手法間の解の一致を検証する DGP**: 処置ユニットの予測子ベクトルがドナー凸包の外に出る設計だと、予測子経路の違いが結果に出ず「一致した」ように見える。因子モデル DGP ＋手計算での突き合わせで判定する
- **OpenMP の並列ループは `num_threads(cs_num_threads())`（`src/threads.h`）で上限を付け、本体を `try/catch` で囲む** — CRAN ポリシーは check 中 2 スレッドまで（`_R_CHECK_LIMIT_CORES_` で自動判定、利用者は `options(coresynth.threads=)`）。例外が OpenMP 領域を抜けると `std::terminate` で R ごと落ちる。`cs_num_threads()` は R の状態を読むので**並列領域の外（メインスレッド）でのみ**呼ぶ。スレッド内では Rcpp 型（`Nullable` 等）も使わない
- **Formula の項は列名そのものに限る**: `scm_fit()` は `all.vars()` で列名を取り `data` から直接読むので、`log(y)` や `I(t > 1990)` を許すと変換が黙って捨てられる。入口で `is.name()` 検査してエラーにしている
- `gh` CLI はインストール済み・認証済み（アカウント `yo5uke`）

---

## ビルド・テスト

手順は `.claude/skills/build-install/SKILL.md`（Windows の DLL ロック解除・サイレントロールバック回避込み）。R は 4.6.1（`bin\` 直下の実行ファイル）。新しい export 関数・S3 method を追加したら、コミット前に `.claude/skills/pkgdown-check/SKILL.md` で `_pkgdown.yml` の reference index を検証する（漏れると pkgdown デプロイが全体停止する）。

---

## コミット

明示的な指示があるまでコミット・プッシュは行わない。作業後は変更を残して報告のみ。

- 件名: Conventional Commits・英語・命令形（`type: summary` / `type(scope): summary`）
- 本文: 英語。「何をしたか」ではなく「なぜ」を書く。自明な変更は本文省略可
- AI ツールの署名・Co-Authored-By は入れない
- バージョン番号はエージェントが判断してよい（破壊的変更・CRAN 提出時期はユーザーに確認）。`Version:` を触る前に **CRAN 現行版を実際に参照**。詳細は `.claude/rules/versioning.md`
- CRAN 提出用 `.tar.gz` の作成・アップロードはユーザーが行う

---

## コーディング規約

- C++ 関数名: `snake_case_cpp` サフィックス。R ラッパー: `fit_<method>_cpp()` が内部関数、`scm_fit()` が公開 API
- 行列の次元規約: `T × N`（行=時間、列=ユニット）統一
- テストは `tests/testthat/test-integration.R` に集約（オブジェクトモデルのみ `test-object-model.R`）
- ggplot の `aes()` 変数は `utils::globalVariables()` に登録する
- コメント: 論文セクション参照は `§` でなく `S.`。開発履歴メモ（いつ・なぜ変更したか）や知らなくても困らないことは書かない。履歴はコミットメッセージ・NEWS.md へ
- `Config/roxygen2/version`（DESCRIPTION）を使う。`RoxygenNote:` が `devtools::document()` で消えるのは roxygen2 8.0.0 以降の正しい移行動作 — **復元しない**

---

## 既知の制限・スコープ外

- **SI**: 共変量未対応（sharp/staggered とも）。原著 Appendix E の時不変共変量の拡張は未実装
- **SCM staggered**: `predictors = NULL` 限定（`covariates` partial-out のみ）
- **plot()**: staggered fit は未対応（明示エラー）。cohort 別の系列は `augment()` の long df を使う
- **MC/TASC の staggered**: `control_group` を選べない。推論関数は sharp/staggered とも無し
- 意図的に対応しない: multivariate outcome の同時推定、SI の time-varying covariates、Callaway & Sant'Anna 型 group-time CATT、Hirshberg & Klosin 型 conformal CI

保留中の項目・未着手の手法は `.claude/docs/backlog.md` 参照。
