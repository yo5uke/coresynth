---
paths:
  - "R/utils.R"
  - "R/scm_fit.R"
---

# 入力検証の設計原則と罠

- **`scm_fit()` で `...` を推定量固有 formal より前に置く**: 後ろに置くと R の部分マッチで MC の `lambda` が SCM の `lambda_pen` に吸われ無言で無効になる。`...` より前の引数は `.check_method_args()` で非対応手法への混入を拒否する（`...` 以降の引数は完全名でしか一致しない）
- **型チェック**: `scm_fit()` は factor/character の outcome・treatment を拒否する（`as.numeric(factor)` はレベルコードを返しデータを無言破壊する）。非整数 `d` も `as.integer()` 切り捨てを防ぐため拒否し、`as.integer(round())` で coerce する
- **sharp 限定の対照ゼロ拒否**: `panel_to_matrices()` がコントロール 0 をエラーにするのは sharp のみ。staggered では never-treated 0 でも MC/TASC が未処置セルから補完できるため通す。再フィット型（SCM/GSC/SI: 2 未満、SDID: 0）は `.check_control_group()` で手法名つきエラー（兼 `"clean"` の非推奨警告）を出す
- **欠損値の扱い**: `.check_panel_complete()`（非有限値エラー）は SCM/SDID/GSC/SI の入口でのみ呼び、MC/TASC では呼ばない（MC は観測マスク `O`、TASC は列平均初期化＋観測期間のみ EM）。一方、`never-treated = 0` 等で識別不能になった期間・ユニットは `.drop_unidentified()` が警告つきで落としてパネルを組み直す（TASC は状態方程式で外挿できるためユニットのみ落とす）
- **処置解除（`1 -> 0`）**: 警告を出したうえで吸収的処置として扱う（データ上の存在セルのみで判定し、採用後の行欠落は解除とみなさない）
- **`panel_to_tensor()` の採用時点判定**: `which(D[, j] > 0L)` を使う（`== 1` だと multi-arm の `d >= 2` が `NA` になり `is_sharp` を誤判定する）
- **`.check_pred_windows()` の役割分担**: パネル外の時点を含む窓（一部重なりでラベルと中身が乖離）と処置後に届く窓には警告を出す。重なりゼロの窓は下流の `build_predictor_matrices()` の非有限値エラーに任せ、二重通知を避ける
