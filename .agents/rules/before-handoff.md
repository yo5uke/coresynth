---
trigger: glob
globs: R/*.R,NEWS.md
---

# 作業を渡す前に振り返ること

- **Rd に載らない引数**: export 関数に `...` を足したときや、`@rdname` を共有するトピックに引数つきのメソッドを足したときは `@param` が抜けやすく、`R CMD check` の WARNING になる（`tools::checkDocFiles(dir = ".")` ですぐ分かる）
- **@examples の総称関数**: `tidy()` などは broom を attach しない check 環境では見つからない。`broom::tidy()` と書く
- **roxygen の数式**: バッククォート内の `\hat\sigma` は文字のまま表示される。数式は `\eqn{}`
- **NEWS.md**: roxygen の `[fn()]` リンク記法は効かない。既存エントリの書き方に合わせる
- **ID の並び順**: unit・cluster の ID は数値・文字列・factor のどれでも来る。`as.character()` してから `sort()` すると `"10" < "2"` になる
- **fit オブジェクトに足すデータ**: 一部の経路でしか使わない行列は、その条件のときだけ持たせる
- **ドキュメントの追従**: export や推論の仕様を足したら、`AGENTS.md` の対応表と該当する `.claude/rules/*.md` に反映が要るか考える
