---
trigger: glob
globs: "R/*.R,src/*.cpp,NEWS.md,.agents/**,.claude/**"
---

# 作業を渡す前に振り返ること

- **速度回帰チェック**: 推定・最適化・推論などの計算経路を更新した後は、関連する既存経路を複数条件（単一フィット、推論/プラセボ、予測子有無、パネル規模など）で実測し、精度とあわせて速度が低下していないことを確認する（`.agents/docs/performance.md` 参照）
- **Rd に載らない引数**: export 関数に `...` を足したときや、`@rdname` を共有するトピックに引数つきメソッドを足したときは `@param` 漏れで `R CMD check` WARNING になる（`tools::checkDocFiles(dir = ".")` で即確認）
- **`@examples` の総称関数**: `tidy()` 等は `broom` を attach しない check 環境では見つからないため、`broom::tidy()` のように名前空間付きで書く
- **roxygen の数式**: バッククォート内の `\hat\sigma` 等は文字のまま展開される。数式は `\eqn{}` を使う
- **NEWS.md**: roxygen の `[fn()]` リンク記法は効かない。既存エントリの記法に合わせる
- **ID の並び順**: unit / cluster の ID は数値・文字列・factor いずれも来うる。`as.character()` 後に `sort()` すると `"10" < "2"` になる点に注意
- **fit オブジェクトの補助データ**: 一部の経路（例: 共変量経路の `Y_raw`/`X_arr`）でのみ使う行列は、該当条件のときだけ保持し、不要経路では `NULL` にする
- **AI ドキュメント同期**: export 関数や推論仕様を変更したら `AGENTS.md` の対応表と該当ルールを更新する。`.agents/rules/`・`.claude/rules/`・`skills/` のいずれかを編集したら、必ずもう一方（フロントマター以外）も同期する
