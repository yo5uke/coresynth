---
name: pkgdown-check
description: 新しいエクスポート関数やS3メソッド（plot.*、print.*等）を追加した後、コミット前にpkgdownデプロイが壊れないか検証する。_pkgdown.ymlのreference indexに載っているかを確認する。
---

新しいエクスポート関数・S3 method を追加したら、コミット前に以下を確認する。
**これを怠ると GitHub Actions の pkgdown デプロイが失敗する**
（`_pkgdown.yml` の reference index にトピックが無いと `build_reference_index()`
がエラーで停止し、サイト全体が更新されない）。

## 手順

1. `devtools::document()` を実行し、`man/*.Rd` が生成されていることを確認する

   ```powershell
   & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "devtools::document()"
   ```

2. **`_pkgdown.yml` の `reference:` セクションに新しいトピックが追加されているか
   確認する**（S3 method は `generic.class` 名で追加。generic と class の対応を
   間違えない）。追加されていなければ `_pkgdown.yml` を編集する

3. ローカルで実際にビルドを確認する（`devtools::document()` も R CMD check も
   これを検出しない）:

   ```powershell
   & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "pkgdown::build_reference_index(pkgdown::as_pkgdown('.'))"
   ```

   エラーなく `Writing `reference/index.html`` まで到達すれば OK。`docs/` は
   gitignore 管理なので生成物はコミット対象外。

4. 逆に、`_pkgdown.yml` に載っているのに `man/` から消えた関数がないかも確認する
   （関数の削除・rename 時に発生）

このチェックは R のみの変更でも必要（`plot.scm_placebo` の追加漏れで pkgdown
デプロイが失敗した実例あり）。
