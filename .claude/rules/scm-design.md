---
paths:
  - "R/scm_design.R"
---

# SCM-Design（Abadie & Zhao 2026）

実験設計段階で処置候補を選ぶための SCM。3 バリアント: `"base"` / `"weakly_targeted"` / `"unit_level"`。推論はブランク期間の順列検定 + スプリット共形 CI。

**eq.(9)/(10) の解法**: 処置候補集合 S を総当たり列挙し、S ごとに (w,v) を最適化する（原著付属コード [SCDesign](https://github.com/jinglongzhao2/SCDesign) と同じ戦略で、原著の Gurobi 非凸 QCQP は使っていない）。`weakly_targeted`（eq.9）は S 固定下で w,v が結合した jointly convex QP になるため **2 ブロック交互最小化**で解く（w-step は標準 SC フィットに帰着する）。`unit_level`（eq.10）は v が w と分離可能なので、per-unit の loss を先に計算してから w の QP に線形項として折り込む（`G = X_S'X_S`、`c = X_S'X_bar − (ξ/2)·loss` の simplex QP を `solve_simplex_qp()` で直接解く。以前の「線形項を目標ベクトルに畳み込む」方式は G が特異だと線形項の一部を落としていた）。

**同点の扱い**: 重み 0 のユニットだけが違う候補集合は目的関数が厳密に同点になる（ξ が大きく w が 1 点に集中すると頻発）。`obj < best_obj` の厳密比較だと丸め誤差の符号で勝者が決まり、BLAS 次第で `treated_units` が変わる（2026-11-19、Ubuntu の oldrel/devel だけでテストが落ちた）。相対 1e-10 未満の差は同点として**列挙順で先の集合**（m の小さい方が先）を残す。`treated_units` は候補集合 S であって w のサポートではない — 「どのユニットに効いているか」を見るテストは `w > 0` で比べること。

**共形 CI の分位点は type 1**（1.1.0 で修正）: 式 (18) は経験分布の下限 `inf{z : F_n(z) ≥ 1−α}`、つまり順序統計量 `|u|_(⌈(1−α)n⌉)` で、`quantile(type = 1)` に当たる。既定の type 7 はその下を補間するので、ブランク期間が少ないと CI が狭くなっていた。順列のサンプリング（組み合わせが 1 万を超えるとき）は固定 seed 42 のまま、呼び出し側の `.Random.seed` を保存・復元する（以前は上書きしていた）。p 値は不変
