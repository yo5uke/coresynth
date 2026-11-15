# 保留中の項目・文献バックログ

着手していない拡張とその理由。着手前に原論文の精読が必須（ここに書いた式・引用は二次情報を含む）。

## 実装が既にあるが未解決の欠落

| 項目 | 内容 |
|---|---|
| not-yet-treated ドナー（`control_group` の第 3 の値） | BFR の `D = {T_i > T(g) + K}` は推定窓 k ≤ K の打ち切りとセット。coresynth は事後窓を末尾まで使うので、対応するには全推論関数が前提にしている「事後行 = T_pre+1..TT」を崩す必要があり、波及範囲が広い |
| TASC の `r` 選択・予測区間 | 論文に d の選択規則は無い（ユーザー指定）。論文は状態共分散から予測区間を出している。`tasc.md` 参照 |
| SI の共変量 | Agarwal et al. Appendix E（時不変共変量を事前アウトカムの下に積む）。`si.md` 参照 |
| sharp SCM の no-op 引数 | `.claude/rules/scm.md` 参照。エラーにすると `control_group` を手法共通で渡しているコードが止まる |

## 利用者のコードや戻り値の意味に触れる変更（ユーザー判断待ち）

`versioning.md` の基準でメジャーに当たるか、当たりうるもの。数値が変わるだけの是正（論文・ドキュメントどおりに近づける修正）はここに積まず実装する。

| 項目 | 現状と判断が要る点 |
|---|---|
| `summary.coresynth()` | summary クラスを返さず表示だけする（CRAN 上は問題なし、慣習上の改善点）。戻り値の構造が変わるので `versioning.md` の基準ではメジャー |
| `scm_fit()` の無名の追加引数 | `.check_method_args()` は `...` の名前なし引数を拒否する文言を持つが、**全部が無名だと `names(list(...))` が NULL になり素通りする**。`scm_fit(f, data, "mc", 5)` の 5 は位置で `lambda` に、`"sdid"` なら `zeta2` に、`"tasc"`・`"gsc"` なら `r` に入る。検査の穴だが、塞ぐと「受理していた入力の拒否」になる |
| 片側 conformal の開いた側 | 対立仮説側は構造的に非有界だが、`ci_upper`（`"greater"`）/`ci_lower`（`"less"`）は初期グリッドの端の値を返している。`Inf`/`-Inf` が正直だが戻り値の意味が変わる |
| SDID の重み QP の条件数 | 理論上はユニット水準に不変だが、時間重みの QP が期間ごとの横断平均しか引かないため、水準が大きいと FISTA の許容誤差の範囲で結果が 4 桁目で動く（`sdid.md`）。ユニットごとの平均も引けば不変になるが、本体の推定値が丸め以上に動く |

### 効率化の余地（結果は不変、未着手）

- **SI の jackknife**: 処置ユニットを落とす回でも対照の SVD をやり直している（GSC は 2026-11-15 に対応済み）。`si_pcr_cpp()` が SVD を内部に持つので、使い回すには C++ 側の分割が要る
- **`sdid_placebo_cpp()`**: `sdid_inference()` が使わなくなった export 関数。次のメジャーで削除候補

## 未着手の手法・機能

- **Cattaneo, Feng, Palomba & Titiunik (2025, REStat)** — staggered SC の非漸近予測区間（conic 最適化）。CRAN に companion の `scpi` が既にあるため、再実装の価値は「6 手法統一インターフェースで完結させる」以外に薄い。現状 staggered で `conformal_inference()` が使えない穴（代わりに `scm_inference()` の wild bootstrap）に当たる
- **Distributional SC (Gunsilius 2023, Econometrica)** — 単位×時間×分位点のデータ構造が要り、`panel_to_matrices()` も Formula も拡張が必要な「7 つ目の柱」。分布データを持つユーザーの需要が確認できるまで保留
- 優先度低: Brown & Butts (2025) の GSC 向け GMM 推論、Li, Shen & Zhou (2024) の factor 型 CI
- Watch（行動不要）: Arkhangelsky & Samkov の Sequential SDID、Lei & Sudijono の refined placebo test、Imbens & Viviano / Arkhangelsky & Hirshberg の confounding 理論
