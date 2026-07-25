---
paths:
  - "R/plot.R"
---

# plot 内部規約

**既知の制限**: staggered fit は `Y_synth = NULL` のため `plot()` 未対応（明示エラー）。cohort 別の系列は `augment()` の long df を使う。

- **`.line_style(default, override)`** は `list` なら `modifyList` マージ、`NULL`/`FALSE` なら NULL（＝非表示）— `geom_vline`/`geom_hline` は一度追加すると除去できないため、抑制判定は**レイヤー追加前**に行う。`colors`/`labels` のユーザー向けキーは 1 単語識別子（`treated`/`synthetic`/`donors`/`placebo`）— `key_map`（識別子→表示名）でマージ後に表示名キーへ変換して ggplot スケールに渡す（キー＝表示文字列だと `labels` での relabel と系列指定が混線するため識別子に分離している）。`.merge_named_vec()` は未知の系列名をエラー（有効名一覧を提示）。trend の `labels` は color/linetype **両スケール**に同一値を渡す（片方だけだと統合凡例が分裂）。plot 引数はすべてデフォルトで既存の見た目を完全再現（後方互換）
- **plot の `align`**: `.align_offset()` が合成系列のプリ期間レベル差を返し、trend/gap 共通で `Y_synth` に加算。SDID 固有の λ 加重整列は `sdid.md` 参照。他手法は単純プリ平均差（ほぼ 0 の微修正）
- **trend の `show_donors`**: `Donors` 系列は color/linetype 両スケールに追加して凡例統合を維持。`breaks` は donors 表示時のみ明示（非表示時に `breaks=NULL` を渡すと凡例自体が消えるので `waiver()`）。ドナー選択は `unit_weights` 降順 → GSC/MC/TASC は明示エラー
- **SDID weights 2 パネル**: `coord_flip` は `facet_wrap(scales="free")` と併用不可のため SDID パスのみ native 横棒で実装。詳細は `sdid.md` 参照
- **`vline_offset`（trend/gap/scm_placebo gaps 共通）**: `.vline_position()` が `T_pre + 1 + offset` を実数インデックスとして扱い、整数部で `times` を引いて小数部を隣接時点差で線形補間 — Date/POSIXct は `times` が既に数値化されているため同じ式で動く。範囲外は警告して線を描かない（`NA` を返し `anyNA(treat_time)` でレイヤーごと抑制）。`vline = list(xintercept = ...)` は `.vline_split()` が style リストから分離し（文字列は軸クラスに coerce）、`vline_offset` との併用はエラー — 位置指定が二重にならないようにするための排他制御
