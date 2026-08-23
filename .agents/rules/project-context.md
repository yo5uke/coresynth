# プロジェクトコンテキスト（coresynth）

このファイルは `CLAUDE.md` の要約。詳細・最新情報は必ずリポジトリルートの `CLAUDE.md` を確認すること（このファイルより `CLAUDE.md` が優先）。

## 概要
**coresynth** は合成コントロール法（SCM）とその派生手法を RcppArmadillo で高速実装した R パッケージ。全推定手法に統一的な Formula インターフェースを提供し、最適化（QP・SVD・カルマンフィルタ）は C++ 実装。

## 対象手法と理論的根拠の参照先
| 手法 | 詳細ファイル |
|---|---|
| SCM | `.claude/rules/scm.md` |
| SDID | `.claude/rules/sdid.md` |
| GSC | `.claude/rules/gsc.md` |
| MC | `.claude/rules/mc.md` |
| TASC | `.claude/rules/tasc.md` |
| SI | `.claude/rules/si.md` |
| SCM-Design | `.claude/rules/scm-design.md` |

手法固有のコードをレビューする際は、必ず対応するファイルを読んでから当たること（アルゴリズム決定の根拠・論文照合がそこにまとまっている）。

横断的ルール:
- `.claude/rules/object-model.md` — fit オブジェクトのクラス構造・アクセサ・staggered 推論の共通ロジック
- `.claude/rules/plot.md` — plot 内部規約
- `.claude/rules/input-validation.md` — 入力検証
- `.claude/rules/versioning.md` — バージョン番号・CRAN 提出方針
- `.claude/rules/conformal.md` — Conformal Inference（`R/conformal.R`、scm/sdid/gsc/mc/si の sharp fit に対応）

随時参照:
- `.claude/docs/theory-crosscheck.md` — 全手法横断の論文照合表
- `.claude/docs/performance.md` — ベンチマーク

## Staggered Adoption 共通仕様（Clarke et al. 2023）
コホート g（採用タイミングが同じ処置ユニット群）ごとに推定した推定量を `N_tr,g * T_post,g` で重み付け平均する。コントロール群は `"clean"`（デフォルト、never-treated + 将来採用者）と `"never_treated"`（統制ユニットのみ）の 2 種。返却構造は `staggered = TRUE`、`cohort_estimates`（data.frame）、`cohort_fits`（list）。数式・詳細は `CLAUDE.md` 参照。

## ファイル構成（要点）
```
coresynth/
├── R/            scm_fit.R(dispatch)・scm.R/sdid.R/gsc.R/mc.R/tasc.R/si.R/scm_design.R(各手法)・
│                 conformal.R・utils.R・accessors.R・broom.R・plot.R・export.R
├── src/          optim.cpp・scm.cpp・sdid.cpp・gsc.cpp・mc.cpp・tasc.cpp・si.cpp・inference.cpp
├── tests/testthat/  test-integration.R(大半のテスト)・test-object-model.R(クラス/アクセサ/dispatch)
├── .claude/      rules/, docs/, skills/ — Claude Code 向け詳細ドキュメント
└── sandbox/, papers/, plan/  gitignore 対象（ベンチマーク・論文MD・計画メモ）
```

## 推論関数の対応表
| 手法 | sharp 推論 | staggered 推論 | `conformal_inference()` |
|---|---|---|---|
| SCM | `mspe_ratio_pval()`（順列検定） | `scm_inference()`（wild bootstrap） | 対応 |
| SDID | `sdid_inference()`（4 method 全て） | 同左 | 対応 |
| GSC | `gsc_boot()` / `gsc_inference()` | `gsc_inference()` | 対応 |
| MC | — | — | 対応 |
| TASC | — | — | 非対応（推論ツール自体なし） |
| SI | `si_inference()`（multi-arm 含む） | 同左 | 対応 |

`jackknife_global` は staggered 専用（sharp・sharp multi-arm では不可）。`conformal_inference()` は sharp fit 限定（staggered/multi-arm/tasc はエラーで手法別関数へ誘導）。
