# 引き継ぎメモ: CRAN 提出前レビュー（2026-10-10, クラウドセッション）

> **一時ファイル**。ローカルで読み込んで作業を再開したら削除してよい（main へのマージ前に消す想定）。
> ブランチ: `ccr-d686b342-k9okjz`（基点 `03e6015`）。`git fetch origin && git checkout ccr-d686b342-k9okjz` で取得。

## 0. 状況サマリ

- 依頼: 1.0.0 を CRAN 提出する前に、潜在バグの検証と、精度を保ったままの高速化・可読性改善のリファクタ
- 実施済み: バグ修正 5 件・高速化 2 件・整理 3 件 → 6 コミットで push 済み（下表）
- **検証は Linux + R 4.3.3（apt）で実施**。R 4.6.1 での検証は**未実施**（クラウド環境では CRAN/Posit が遮断されていて、conda-forge にも R 4.6 向けの Rcpp 等が無く中断）→ **ローカル（Windows + R 4.6.1）での再検証が最優先**
- `Version:` は**未変更**（CRAN へ到達できず、現行版を `versioning.md` の手順どおりに参照できなかったため）。NEWS は既存の `# coresynth 1.0.0` 見出しに追記済み

| コミット | 内容 |
|---|---|
| `6d986f5` | fix: `gsc_boot()` の CI を ATT 中心に修正 / Formula の変換項（`log(y)` 等）をエラー化 |
| `ad9a0bf` | fix: OpenMP のスレッド上限（CRAN ポリシー）と並列ループ内の例外捕捉 |
| `ed6689e` | perf: SDID の重み QP を因子形式で解く / TASC のカルマンゲインを r×r で解く / 未使用の `solve_simplex_qp_lr()` 削除 |
| `2993d29` | refactor: GSC EM の不変行列をループ外へ / SI-PCR のループを行列積に |
| `643d5b3` | test: `tests/testthat.R` 追加＋回帰テスト 12 件 |
| `16870c8` | docs: NEWS・CLAUDE.md・rules（gsc/sdid/tasc）に記録 |

## 1. ローカルでやること（チェックリスト）

1. [ ] `devtools::document()` → **差分が出ないこと**を確認（`man/gsc_boot.Rd` と `man/coresynth-package.Rd` はクラウドに roxygen2 8.1.0 が無かったため**手編集**。roxygen 出力と書式がずれていれば、生成された方を採用する）
2. [ ] `Rcpp::compileAttributes()` → 差分なし（クラウドでは確認済み）
3. [ ] `build-install` skill（C++ を変更しているので `.o` を消してから INSTALL、`DONE (coresynth)` と DLL の更新を確認）
4. [ ] `devtools::test()` → 期待値は **1883 pass / 0 fail / 2 skip**（R 4.3.3 での結果。修正前のベースラインは 1871 pass）
5. [ ] `devtools::check(args = c("--no-manual", "--run-donttest", "--as-cran"))` → 0/0/0 を確認。**今回から check 内でテストが走る**（約 60 秒、CPU/経過 ≈ 1.05）
6. [ ] 速度の確認: `sandbox/24_bench_vs_synth_large.R`（SCM 側。今回 SCM の内側 QP は触っていないので桁が変わらないこと）＋付録 A のベンチ
7. [ ] 付録 A のベンチで、推定値が「修正前」の列と 10 桁一致すること（MC だけは元々非決定ではないので一致するはず。下表参照）
8. [ ] バージョン番号を決める（§5）→ `chore(release)` コミット
9. [ ] このファイルを削除

## 2. 修正したバグ（詳細）

### 2.1 `gsc_boot()` の CI が 0 中心だった（`R/gsc.R`）
- 複製は効果ゼロで生成しているので、ATT* は推定誤差 δ̂−δ の分布（Xu 2017 Algorithm 2）
- 旧: `ci = quantile(boot, c(α/2, 1−α/2))` → 0 中心の**受容域**。再現データでは ATT=1.82 に対して CI=[−0.28, 0.30]
- 新: `ci_lower = ATT − q(1−α/2)`, `ci_upper = ATT − q(α/2)`。`p_value`・`se`・`boot_dist` は不変
- README の「95% CI」表示と `export_json()` に影響。`vignettes/articles/inference.Rmd` に 1 文追記、`man/gsc_boot.Rd` を手編集
- 規約は `.claude/rules/gsc.md` に記録済み

### 2.2 Formula の変換が黙って捨てられていた（`R/scm_fit.R`）
- `all.vars()` で列名を取って `data` から直接読むため、`log(y) ~ D | id + t` が `y` のまま推定されていた（推定値が完全に同一）。`y ~ I(t > T0) | ...` では `t` 列がそのまま処置変数になっていた
- 入口で LHS と処置項に `is.name()` 検査を入れ、変換項はエラーにした（`formula(f_parts, lhs=1, rhs=0)[[2L]]` / `formula(f_parts, lhs=0, rhs=1)[[2L]]`）
- 単位・時間の項（`| id + time`）は従来どおり（変換されても列名を読むだけで実害がないため）

### 2.3 OpenMP のスレッド上限（`src/threads.h` 新規、`src/inference.cpp`、`src/scm.cpp`）
- CRAN ポリシーでは check 中は 2 スレッドまで。旧実装は全コアを使っていた
- `cs_num_threads()` の優先順位: `options(coresynth.threads = n)`（正の数）→ `_R_CHECK_LIMIT_CORES_` が設定されていれば（`--as-cran`）2 → OpenMP の既定値（`OMP_NUM_THREADS`/`OMP_THREAD_LIMIT` を尊重）
- 全 7 箇所の `#pragma omp parallel for` に `num_threads(n_threads)` を付与。`scm_multistart_core` は `spawn_team` のときだけ `cs_num_threads()` を呼ぶ（**R の状態を読むのでメインスレッド限定**）
- 実測の CPU/経過比: 既定 2.93（4 コア）/ `--as-cran` 時 1.98 / option=1 で 0.99
- 結果がスレッド数に依存しないことはテストで確認（`expect_identical`）
- `R/coresynth-package.R` に `@section Parallel computation:` を追加し、Rd を手編集

### 2.4 並列ループ内の例外で R が落ちる（`src/inference.cpp`）
- `sdid_placebo_cpp` と `scm_placebo_cpp`（outcomes-only）に `try/catch` が無かった。例外が OpenMP 領域を抜けると `std::terminate` で R セッションごと落ちる
- 失敗した placebo は NaN にする（`scm_placebo_x_cpp` の既存の契約に合わせた）

### 2.5 `tests/testthat.R` が存在しなかった
- git 履歴上一度も存在していない → **ローカルの check でも GitHub Actions でもテストが走っていなかった**
- 標準形で追加。テストは check 内で約 60 秒

## 3. 高速化（精度維持を確認済み）

| 処理 | 修正前 | 修正後 | 一致 |
|---|---|---|---|
| SDID placebo（N=300） | 4.84s | 0.78s | 推定・SE・p 値が 10 桁一致 |
| SDID jackknife（N=300） | 6.38s | 0.81s | 同上 |
| SDID bootstrap 50 回（N=300） | 1.06s | 0.20s | 同上 |
| SDID staggered jackknife_global（N=150） | 2.01s | 0.83s | 同上 |
| TASC（N=200） | 3.93s | 0.34s | 同上 |
| TASC staggered（N=120） | 0.92s | 0.24s | 同上 |

（クラウドの 4 コア、R 4.3.3。修正前の計測時は裏で check が走っていたので絶対値は参考程度。同じ回で交互に測り直すこと）

- **SDID**（`src/optim.cpp`、`src/sdid.cpp`）: `solve_simplex_ridge_ls(A, b, ridge)` を新設。`min ½‖Aw‖² + ½ρ‖w‖² − (A'b)'w` を、A が横長なら N×N の Gram を作らずに勾配 `A'(Ay) + ρy − c`、Lipschitz 定数 `λ_max(AA') + ρ` で FISTA。反復列は密形式と丸め誤差を除いて同一。FISTA 本体は `fista_simplex_core<Grad>` テンプレートに一本化し、`solve_simplex_qp()` は従来と同じ演算（`Q*y - c`）のまま。新ソルバは Rcpp 型を使わないのでスレッド安全
- **TASC**（`src/tasc.cpp`）: R が対角かつ正なら push-through 恒等式 `K = (I + P W'R⁻¹W)⁻¹ P W'R⁻¹` で r×r の solve にし、`KRK'` は列スケーリングで計算。一般の R や solve 失敗時は従来の直接形にフォールバック。テストでは非対角成分に 1e-300 を入れて直接形を強制し、比較している
- **整理（速度は不変、結果も同一）**: GSC EM の `X'X` をループ外へ出し、`X'R` は `accu(slice % R)` で計算（スタック行列を廃止）。SI-PCR は `W = Pinv * Y_pre_tr`、未使用の `solve_simplex_qp_lr()` は削除（0.4.1 の因子形式化以降どこからも呼ばれていなかった）

## 4. 試して取り消したもの（再挑戦しないこと）

- **MC の収束判定を L 単独から完成行列 M = L + FE に変更 → 取り消し**
  - 動機: λ が大きく L=0 になると、FE が未収束のまま止まるのではと考えた
  - 実測: L=0 でも FE は exact TWFE と 1e-13 まで一致していて、問題は存在しなかった
  - 逆効果: M はレベル（~10）を含んで相対変化が小さく見えるので**早く止まり、精度が悪化**（tol=1e-12 の解からのずれが 6.9e-5 → 4.7e-4）
- MC の CV のウォームスタート: `mc.md` の「不採用」の判断を尊重して未着手

## 5. バージョン番号の推奨

- NEWS から読む限り CRAN 現行版は 0.5.0（**未確認** — ローカルで https://cran.r-project.org/package=coresynth を参照すること）
- 1.0.0 は CRAN 未公開なので、そのまま出しても合法。ただし GitHub の 1.0.0 は 2 週間公開されていて、今回 `gsc_boot()` の出力が変わり、Formula がエラー化した → `versioning.md` の「同じ番号で違う中身」に該当
- **推奨: 今回の修正込みで、1.0.1 として最初の提出を行う**（1.0.0 は CRAN 欠番でよい）。「1.0.0 をそのまま出してから 1.0.1」は、既知の誤った CI を含む版を出すことになり、CRAN の更新間隔（1〜2 ヶ月）にも引っかかるので非推奨
- 1.0.1 にするなら NEWS の見出しを分ける（今回の項目を `# coresynth 1.0.1` に移す）か、1.0.0 の見出しをそのまま 1.0.1 に改名するかはユーザー判断

## 6. ユーザー判断待ち（コード未変更・推定値や推論結果が変わる設計選択）

1. **SDID placebo**: λ と ζ を元のフィットの値で固定し、N_tr>1 でも疑似処置を 1 ユニットにしている。Clarke et al. (2024) Algorithm 4 と synthdid は λ を再推定し、N_tr 個のユニットを疑似処置にする → 現状は N_tr>1 で SE が保守的（大きめ）。`sdid.md` の「λ・ζ² 固定」は意図的な記述あり
2. **SDID bootstrap**: 対照群だけをリサンプル（処置は固定）。Clarke Algorithm 2 は全ユニットをリサンプルする
3. **`gsc_boot()` の p 値**: `mean(|boot| >= |att|)` で +1 補正が無く、ちょうど 0 になりうる（Davison–Hinkley の `(1+#)/(B+1)` が標準）
4. **処置の途中解除**（d が 1→0 に戻る）: 吸収的とみなして黙って処理される（6 手法とも推定値は解除なしと同一）。`panel_to_matrices()` で警告またはエラーにする余地あり
5. `summary.coresynth()` は summary クラスを返さずに表示だけする（CRAN 的には問題なし、慣習上の改善点）

## 7. 確認したが問題なかったもの

- `R CMD check --as-cran`（修正前・修正後とも）: 残った WARNING と NOTE は環境由来のみ（ネット遮断で CRAN 索引に到達不可、`en_US.UTF-8` ロケール無し、Ubuntu の R が付ける `-mno-omit-leaf-frame-pointer`、デバッグシンボルで libs が 11.5MB）
- 境界ケース（6 手法すべて、推定値が基準と一致）: id が factor / 時間が Date / T_post=1 / 処置ユニットの id が末尾 / 行シャッフル / logical の D / tbl_df クラス
- R コード中の非 ASCII はコメントのみ（plot タイトルの規約 OK）。`paste0("Donor ", ...)` も無し
- `1:length()` 系・`sample()` の長さ 1 の罠・`T`/`F` を値として使う箇所: 実害なし（`gsc_boot()` 内のローカル変数 `T <- nrow(F_hat)` は TRUE を隠すが参照はなく無害。気になればリネーム）
- Armadillo の警告がスレッドから出ていないこと（テストログに 0 件）

## 付録 A: ベンチスクリプト

修正前のインストールと修正後のインストールを `R_LIBS` で切り替え、同じ回で交互に走らせる。`Rscript bench.R out.rds`。

```r
suppressMessages(library(coresynth))
args <- commandArgs(TRUE); out <- args[1]
options(coresynth.threads = 4L)
mk <- function(N, TT, T0, n_tr = 1, stag = FALSE, seed = 1) {
  set.seed(seed)
  F <- cbind(1, cumsum(rnorm(TT)), sin(seq_len(TT) / 3))
  L <- cbind(rnorm(N, 3), rnorm(N), rnorm(N))
  Y <- F %*% t(L) + matrix(rnorm(TT * N, sd = .5), TT, N) + 10
  df <- data.frame(id = rep(1:N, each = TT), t = rep(1:TT, N), y = as.vector(Y))
  adopt <- rep(Inf, N)
  adopt[seq_len(n_tr)] <- if (stag) rep(c(T0, T0 + 4, T0 + 8), length.out = n_tr) else T0 + 1
  df$D <- as.integer(df$t >= adopt[df$id])
  df$x <- rnorm(nrow(df)) + df$t / 10
  df$y <- df$y + 2 * df$D + 0.5 * df$x
  df
}
tm <- function(expr) { g <- gc(); t <- system.time(v <- force(expr))[["elapsed"]]; list(t = t, v = v) }
res <- list()
add <- function(name, r, val) { res[[name]] <<- list(t = r$t, val = val); cat(sprintf("%-28s %7.2fs  %s\n", name, r$t, paste(signif(val, 10), collapse = " "))) }

d1 <- mk(300, 40, 30, n_tr = 3)
r <- tm(scm_fit(y ~ D | id + t, d1, method = "sdid")); add("sdid_fit N300", r, r$v$estimate)
f_sdid <- r$v
r <- tm(sdid_inference(f_sdid, "placebo")); add("sdid_placebo N300", r, c(r$v$se, r$v$p_value))
r <- tm(sdid_inference(f_sdid, "jackknife")); add("sdid_jack N300", r, r$v$se)
r <- tm(sdid_inference(f_sdid, "bootstrap", n_boot = 50, seed = 1)); add("sdid_boot50 N300", r, r$v$se)
d2 <- mk(150, 40, 25, n_tr = 6, stag = TRUE)
r <- tm(scm_fit(y ~ D | id + t, d2, method = "sdid")); add("sdid_stag N150", r, r$v$estimate)
r <- tm(sdid_inference(r$v, "jackknife_global")); add("sdid_stag_jkg N150", r, r$v$se)

d3 <- mk(200, 50, 35, n_tr = 2)
r <- tm(scm_fit(y ~ D | id + t, d3, method = "tasc")); add("tasc N200", r, r$v$estimate)
d3s <- mk(120, 50, 30, n_tr = 6, stag = TRUE)
r <- tm(scm_fit(y ~ D | id + t, d3s, method = "tasc")); add("tasc_stag N120", r, r$v$estimate)
r <- tm(scm_fit(y ~ D | id + t, d3, method = "mc")); add("mc N200", r, c(r$v$estimate, r$v$lambda))
r <- tm(scm_fit(y ~ D | id + t, d3s, method = "mc")); add("mc_stag N120", r, c(r$v$estimate, r$v$lambda))
r <- tm(scm_fit(y ~ D | id + t, d3, method = "gsc")); add("gsc N200", r, c(r$v$estimate, r$v$r))
r <- tm(scm_fit(y ~ D | id + t, d3, method = "gsc", covariates = "x")); add("gsc_cov N200", r, c(r$v$estimate, r$v$beta))
f_gsc <- r$v
r <- tm(gsc_inference(f_gsc, "bootstrap", n_boot = 50, seed = 1)); add("gsc_cov_boot50", r, r$v$se)
r <- tm(scm_fit(y ~ D | id + t, d3, method = "si")); add("si N200", r, r$v$estimate)
d4 <- mk(80, 40, 30)
r <- suppressMessages(tm(scm_fit(y ~ D | id + t, d4, method = "scm"))); add("scm N80", r, r$v$estimate)
r <- tm(mspe_ratio_pval(r$v)); add("scm_placebo N80", r, r$v$p_value)
saveRDS(res, out)
```

修正前の出力（R 4.3.3、4 コア）:

```
sdid_fit N300                   0.06s  2.018887379
sdid_placebo N300               4.84s  0.2800276024 0.003355704698
sdid_jack N300                  6.38s  0.05326316693
sdid_boot50 N300                1.06s  0.05532946361
sdid_stag N150                  0.03s  2.087698546
sdid_stag_jkg N150              2.01s  0.04239027047
tasc N200                       3.93s  2.159555348
tasc_stag N120                  0.92s  1.948157528
mc N200                         2.55s  2.106221997 8.994500938
mc_stag N120                    2.68s  1.974226311 9.886088144
gsc N200                        0.04s  2.081754539 3
gsc_cov N200                    0.55s  2.006845632 0.4972457613
gsc_cov_boot50                  4.61s  0.01518929808
si N200                         0.03s  2.084957163
scm N80                         0.06s  1.640826264
scm_placebo N80                 0.92s  0.0125
```

最終版（MC の変更を取り消した後のビルド）での実測: 全 16 項目で値の最大差は 1.9e-14 以下（MC・GSC 無共変量・SI・SCM は完全一致の 0）。時間は SDID placebo 0.75s・jackknife 0.86s・TASC 0.31s、それ以外は計測誤差の範囲。

## 付録 B: クラウド環境の制約（再発時の参考）

- CRAN（cran/cloud.r-project.org）・Posit（cdn.posit.co, packagemanager）・github.com のリリースはプロキシで遮断。conda-forge と archive.ubuntu.com には到達できる
- conda-forge には r-base 4.6.1 があるが、R パッケージの 4.6 向けビルドはほぼ無い（2026-10 時点）。Ubuntu の `dists/devel/universe/source/Sources.xz` から CRAN のソース tarball（`r-cran-*_*.orig.tar.gz`）を集めて、ローカル CRAN リポジトリを作る案は、ユーザーの指示で中断した
