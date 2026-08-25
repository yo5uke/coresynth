# Git 運用ルール（coresynth）

Claude Code に適用しているルールと同一。

- **明示的な指示があるまでコミット・プッシュは行わない**。作業後は変更を作業ツリーに残したまま、変更内容を報告するだけにとどめる。
- コミットメッセージの書き方:
  - 件名: Conventional Commits 形式・英語・命令形（`type: summary` または `type(scope): summary`）。主な type は `feat`/`fix`/`perf`/`docs`/`chore`（`chore(release)` を含む）/`ci`
  - 本文: 英語。「何をしたか」ではなく「なぜ」（動機・原因・トレードオフ・検証結果）を書く。自明な変更は本文を省略してよい
  - AI ツールの署名・Co-Authored-By は入れない
- force push・`git reset --hard`・`git checkout --`/`restore` による変更破棄など、破壊的・元に戻しにくい操作は明示的な指示なしに行わない。
- 既存コミットへの amend より新規コミットを優先する（明示的に amend を指示された場合を除く）。
- pre-commit hook 等は `--no-verify` でスキップしない。フックが失敗したら原因を調べて修正してから再コミットする。
- 破壊的操作の前（checkout/restore/reset/clean 等）には必ず `git status` で作業ツリーの状態を確認し、未コミットの変更があれば先にユーザーに確認する。
- ステージング・コミット前に `git status`/`git diff` で内容を確認し、意図しないファイル（秘密情報を含むファイル等）が混ざっていないか確認する。
- CRAN 提出用の `.tar.gz` の作成・アップロードはユーザーが自分で行う。エージェントは行わない（`R CMD build` 等でパッケージ化する作業はスコープ外）。
