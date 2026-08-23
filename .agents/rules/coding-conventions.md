# コーディング規約と既知の制限（coresynth）

## 命名・構造
- C++ 関数名: `snake_case_cpp` サフィックス（例: `sdid_unit_weights_cpp`）
- R ラッパー: `fit_<method>_cpp()` が内部関数、`scm_fit()` が公開 API
- 行列の次元規約: `T × N`（行=時間、列=ユニット）で統一
- テストは `tests/testthat/test-integration.R` に集約（オブジェクトモデルのみ `test-object-model.R`）
- ggplot の `aes()` 変数は `utils::globalVariables()` に登録する（R CMD check NOTE 回避）

## コメント規約
- 論文セクション参照は `§` ではなく `S.` を使う
- 開発履歴メモ（いつ・なぜ変更したか、開発フェーズへの言及など）はコメントに書かない。コミットメッセージ・NEWS.md に書く
- 知らなくても困らないことは書かない（蛇足コメントを避ける）

## roxygen2（バージョン 8.0.0）
```
Config/roxygen2/markdown: TRUE
RoxygenNote: 8.0.0
```
`devtools::document()` 後に `RoxygenNote:` が 7.x に書き戻された場合は即座に `8.0.0` に修正する。

## 新規 export 追加時
新しい export 関数・S3 method（`plot.*`、`print.*` 等）を追加した場合は `_pkgdown.yml` の reference index にトピックがあるか確認する（無いと GitHub Actions の pkgdown デプロイがサイト全体停止する）。

## 既知の制限
- SI の共変量: 未対応（sharp/staggered とも）
- SCM staggered: `predictors = NULL` 限定（`covariates` partial-out のみ）
- `plot()`: staggered fit は `Y_synth = NULL` のため未対応

## スコープ外（意図的な制限）
- multivariate outcome（同時推定）
- SI の time-varying covariates
- Callaway & Sant'Anna (2021) 型 group-time CATT
- Hirshberg & Klosin (2023) 型 conformal CI

詳細な判断根拠は `CLAUDE.md` 本体を参照。
