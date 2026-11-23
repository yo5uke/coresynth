---
trigger: glob
globs: "R/plot.R"
---

# plot 内部規約

- **staggered fit は `plot()` 未対応**: `Y_synth = NULL` のため明示エラーとし、`augment()` のコホート別 long df へ誘導する
- **参照線の抑止はレイヤー追加前に判定**: `geom_vline`/`geom_hline` は追加後に除去できないため、`.line_style()` が `NULL`/`FALSE` で `NULL` を返し、追加自体をスキップする
- **1 単語識別子キーと `key_map`**: `colors` / `labels` / `linetypes` のユーザー向けキーは `treated` / `synthetic` / `donors` / `placebo` とし、マージ後に `key_map` で表示名へ変換する（キーと表示文字列を混用すると `labels` のリラベル時に系列指定と衝突するため）。`trend` の `labels` は color と linetype 両スケールに同一値を渡して統合凡例を維持する
- **`linetypes` の数値変換**: マージ先が character ベクトルのため、`.as_linetype()` で数値コード `0:6` を名前に変換してから `.merge_named_vec()` に渡す（生の数値を混ぜると ggplot が hex ダッシュパターンと誤認する）。`scm_placebo` gaps は既定が両系列 solid でも linetype スケールを常に追加して凡例を統合する
- **trend の `show_donors`**: 非表示時に `breaks = NULL` を渡すと凡例自体が消えるため `waiver()` を使う。ドナー上位選択は `unit_weights` 降順を使うため、重みを持たない GSC/MC/TASC では明示エラーにする
- **`plot.scm_placebo(type = "ratios")`**: サブタイトルが `n_placebo_used` を出すため `mspe_prune` を同じように適用し、`statistic = "post_mspe"` では事後 MSPE を描く。`plot_data()` は剪定を適用しつつ `mspe_threshold` 由来の非有限行は残して報告する
- **SDID weights 2 パネル**: `coord_flip()` は `facet_wrap(scales = "free")` と併用できないため SDID のみ native 横棒で描画する。`align` の SDID 固有定義は `sdid` ルール参照
- **`vline_offset` と `.VLINE_OFFSET_DEFAULT = -1`**: `.vline_position()` が `T_pre + 1 + offset` を線形補間する（Date/POSIXct も同じ式）。範囲外警告および `vline = list(xintercept = ...)` との排他エラーは必ず定数 `.VLINE_OFFSET_DEFAULT`（リテラル `0` ではない）と比較して判定する
