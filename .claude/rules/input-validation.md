---
paths:
  - "R/utils.R"
  - "R/scm_fit.R"
---

# 入力検証

- **`panel_to_matrices()`**（fit 全手法の共通入口）: ①0 行データ ②NA の id/time ③NA の処置指標（従来は該当ユニットが無言で脱落）④負の d ⑤処置ユニットなし ⑥sharp 時コントロールなし、を明示エラー。⑥が **sharp 限定**なのが重要 — staggered で全ユニットが最終的に処置されても MC/TASC は未処置セルから補完できる。コホート再フィット型（SCM/SDID/GSC/SI）の never-treated 2 未満は `.check_control_group()` が手法名つきで明示エラーにする（同関数は `"clean"` の非推奨警告も出す）。`scm_design()` にも 0 行チェック
- **`.check_panel_complete()`**: SCM/SDID/GSC/SI の fit 入口で `!is.finite(Y)` セルをユニット・時点つきで明示エラー（従来は `eig_sym(): decomposition failed` 等の C++ エラー）。欠損対応の MC/TASC は呼ばない — MC は O マスクで、TASC は初期化のみ列平均補完（EM 本体は NA のままカルマン + per-unit M-step）で欠損を扱う
- **`scm_fit()` 型検証**: factor/character の outcome・treatment を拒否（`as.numeric(factor)` はレベルコードを返しデータを無言破壊）。非整数 d も拒否（`as.integer` 切り捨てで 0.9→0 の罠）、coerce は `as.integer(round())`
- **`build_predictor_matrices()`**: X0/X1 の非有限値を predictor 名つきで即エラー（`pred()` の時間窓にデータがないケースで下流 QP が壊れるのを防止）
- **`.check_pred_windows()`**（sharp SCM の予測子経路のみ）: ①パネルに無い時点を含む窓（重なりが一部でもあれば集約は成功してしまい、ラベルだけが嘘になる）②処置後期間に届く窓（処置で動いた値を集約する）を警告。重なりゼロの窓は上記の非有限エラーに任せる（二重通知を避ける）
- **`panel_to_tensor()` の T_adopt**: 採用時点判定は `which(D[, j] > 0L)`（`== 1` だと arm≥2 のユニットが NA になり is_sharp 誤判定）
- **`scm_fit()` の引数検証**（`.check_method_args()`）: ①手法非対応の formal（`predictors`/`nu`/`fixedeff` 等は scm 限定、`covariates` は scm/sdid/gsc）を渡したらエラー ②`...` に当該手法の fit 関数の formal に無い名前があればエラー（誤記・`control_group` を mc/tasc に渡す等）。**`...` は推定量固有 formal より前に置く** — 後ろに置くと R の部分マッチで `lambda`（mc）が `lambda_pen`（scm 限定）に吸われ、MC のペナルティが無言で届かなくなる。`...` より後ろの引数は完全名でしか一致しない性質を利用している
