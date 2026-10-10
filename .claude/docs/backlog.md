# 保留中の項目・文献バックログ

着手していない拡張とその理由。着手前に原論文の精読が必須（ここに書いた式・引用は二次情報を含む）。

## 実装が既にあるが未解決の欠落

| 項目 | 内容 |
|---|---|
| not-yet-treated ドナー（`control_group` の第 3 の値） | BFR の `D = {T_i > T(g) + K}` は推定窓 k ≤ K の打ち切りとセット。coresynth は事後窓を末尾まで使うので、対応するには全推論関数が前提にしている「事後行 = T_pre+1..TT」を崩す必要があり、波及範囲が広い |
| GSC の加法固定効果、TASC の `r` 選択 | 各手法の rules ファイル参照 |
| `augment_scm()` の staggered エラー文言、sharp SCM の no-op 引数 | `.claude/rules/scm.md` 参照 |

## 推定値・推論結果が変わる設計選択（ユーザー判断待ち）

どれも推定値か推論結果が変わるので、実装前に方針を確認する。

| 項目 | 現状と論文との差 |
|---|---|
| SDID placebo | λ・ζ を元のフィットの値で固定し、N_tr>1 でも疑似処置を 1 ユニットにしている（`sdid.md` の「λ・ζ² 固定」は意図的）。Clarke et al. (2024) Algorithm 4 は λ を再推定し N_tr 個を疑似処置にする → N_tr>1 で SE が保守的 |
| SDID bootstrap | 対照群だけをリサンプル（処置は固定）。Clarke Algorithm 2 は全ユニットをリサンプル |
| `gsc_boot()` の p 値 | `mean(\|boot\| >= \|att\|)` で +1 補正が無くちょうど 0 になりうる（Davison–Hinkley の `(1+#)/(B+1)` が標準） |
| 処置の途中解除（d が 1→0） | 吸収的とみなして黙って処理される（6 手法とも推定値は解除なしと同一）。`panel_to_matrices()` で警告かエラーにする余地 |
| `summary.coresynth()` | summary クラスを返さず表示だけする（CRAN 上は問題なし、慣習上の改善点） |

## 未着手の手法・機能

- **Cattaneo, Feng, Palomba & Titiunik (2025, REStat)** — staggered SC の非漸近予測区間（conic 最適化）。CRAN に companion の `scpi` が既にあるため、再実装の価値は「6 手法統一インターフェースで完結させる」以外に薄い。現状 staggered で `conformal_inference()` が使えない穴（代わりに `scm_inference()` の wild bootstrap）に当たる
- **Distributional SC (Gunsilius 2023, Econometrica)** — 単位×時間×分位点のデータ構造が要り、`panel_to_matrices()` も Formula も拡張が必要な「7 つ目の柱」。分布データを持つユーザーの需要が確認できるまで保留
- 優先度低: Brown & Butts (2025) の GSC 向け GMM 推論、Li, Shen & Zhou (2024) の factor 型 CI
- Watch（行動不要）: Arkhangelsky & Samkov の Sequential SDID、Lei & Sudijono の refined placebo test、Imbens & Viviano / Arkhangelsky & Hirshberg の confounding 理論
