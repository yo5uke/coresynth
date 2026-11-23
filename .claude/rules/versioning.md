---
paths:
  - "DESCRIPTION"
  - "NEWS.md"
---

# バージョン運用

番号の決定はエージェントに委任されている。ただし **破壊的変更（major）** が絡む場合と **CRAN 提出時期** を動かす判断はユーザーに確認する。

## 2 つのチャンネルと確認手順

| チャンネル | 実体 | 確認先 |
|---|---|---|
| 開発版 | GitHub `main`（`pak::pak("yo5uke/coresynth")`） | `DESCRIPTION` の `Version:` |
| リリース | CRAN | `https://cran.r-project.org/package=coresynth`（推測せず必ず参照） |

`DESCRIPTION` の版は GitHub に push した時点で「開発版として公開済み」になる。`Version:` を触る前に必ず以下を確認する:
1. CRAN 現行版と公開日を参照する
2. `DESCRIPTION` の現在値が CRAN 公開済みかを見る:
   - **未公開（GitHub 先行中）**: 番号は据え置き、`NEWS.md` も同じ見出しに追記する（開発中に都度刻み直さない）
   - **公開済み**: 下記基準で新しい番号を決める
3. CRAN 提出を伴う場合、直近公開日から 1〜2 ヶ月空いているか確認し、空いていなければ提出時期をユーザーに相談する

## 番号決定基準（基準は「利用者のコードが壊れるか」）

- **patch (`x.y.Z`)**: バグ修正・性能改善・内部変更（API と引数は不変）。推定値や推論結果が変わっても、論文・ドキュメントどおりの推定量に近づける是正は patch
- **minor (`x.Y.0`)**: 後方互換な新機能・新引数・新 export、および既定挙動の変更（既定引数や既定推定量の変更）。必ず `.0` から始める（`0.1.2 -> 0.2.1` のように `.0` を飛ばさない）
- **major (`X.0.0`)**: 既存コードが動かなくなる変更のみ（引数・関数の削除や改名、戻り値構造の変更、受理していた正当入力の拒否。データ誤りを指摘するエラー追加は除く）

数値が変わる修正では、どの版でも `NEWS.md` に **変わる条件と変わる量** を明記する（旧挙動を残す引数は原則作らない）。

## `.9000` 開発版サフィックスと CRAN 制約

- 未公開の据え置き中に中身がまとまって変わる、または次の提出番号が現在値では足りないと確定した場合は、直近で GitHub に出した番号に `.9000` を付け、`NEWS.md` 見出しを `# coresynth (development version)` にする。提出時は `.9000` を落として戻すのではなく**次の番号へ上げる**（戻すとダウングレードになる）
- CRAN 提出版は現行版より成分ごとの数値比較で厳密に大きい必要がある（欠番は可、公開済み番号の再利用・上書きは不可、却下で未公開なら同番号再提出可）

## リリースチェックリスト

1. `NEWS.md` の見出しを確定番号にし、内容が実装と一致するか確認
2. `DESCRIPTION` の `Version:` を更新（`Config/roxygen2/version` が現行と一致しているか確認。`RoxygenNote:` は復元しない）
3. `devtools::document()` → `build-install` スキル → `devtools::test()`
4. `devtools::check(args = c("--no-manual", "--run-donttest"))` で 0 errors / 0 warnings / 0 notes
5. 新規 export / S3 method があれば `pkgdown-check` スキル
6. コミットは「機能」と「`chore(release): bump version to x.y.z`（DESCRIPTION + NEWS.md）」に分ける（`.tar.gz` 作成・CRAN 提出はユーザー担当）
