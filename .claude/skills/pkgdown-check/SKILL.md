---
name: pkgdown-check
description: 新しいエクスポート関数やS3メソッド（plot.*、print.*等）を追加した後、コミット前にpkgdownデプロイが壊れないか検証する。_pkgdown.ymlのreference indexに載っているかを確認する。
---

新しいエクスポート関数・S3 method を追加・削除・rename したら、コミット前に以下を確認する。
`_pkgdown.yml` の reference index にトピックが無いと `build_reference_index()` がエラーで停止し、GitHub Actions の pkgdown デプロイ全体が失敗する（`devtools::document()` も `R CMD check` もこれを検出しない）。

## 手順

1. `devtools::document()` を実行し、`man/*.Rd` の生成を確認する:

   ```powershell
   & "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" -e "devtools::document()"
   ```

2. `_pkgdown.yml` の `reference:` セクションに新しいトピックが追加されているか確認する（S3 method は `generic.class` 名で記載）。逆に、`man/` から消えた関数が `_pkgdown.yml` に残っていないかも確認する。

3. ローカルで reference index のビルドを検証する:

   ```powershell
   & "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" -e "pkgdown::build_reference_index(pkgdown::as_pkgdown('.'))"
   ```

   エラーなく `Writing 'reference/index.html'` まで到達すれば OK（`docs/` は gitignore 管理のためコミット対象外）。
