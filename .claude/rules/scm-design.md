---
paths:
  - "R/scm_design.R"
---

# SCM-Design（Abadie & Zhao 2026）

実験設計段階で処置候補を選ぶための SCM。3 バリアント: `"base"` / `"weakly_targeted"` / `"unit_level"`。推論はブランク期間の順列検定 + スプリット共形 CI。

**eq.(9)/(10) の解法**: 処置候補集合 S を総当たり列挙し、S ごとに (w,v) を最適化する（原著付属コード [SCDesign](https://github.com/jinglongzhao2/SCDesign) と同じ戦略で、原著の Gurobi 非凸 QCQP は使っていない）。`weakly_targeted`（eq.9）は S 固定下で w,v が結合した jointly convex QP になるため **2 ブロック交互最小化**で解く（w-step は標準 SC フィットに帰着する）。`unit_level`（eq.10）は v が w と分離可能なので、per-unit の loss を先に計算してから w の QP に線形項として折り込む。
