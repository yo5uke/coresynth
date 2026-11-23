---
trigger: glob
globs: "R/tasc.R,src/tasc.cpp"
---

# TASC（Rho et al. 2026）

$$z_{t+1} = A z_t + C + \eta_t, \quad y_t = \mu + W z_t + \varepsilon_t$$

- **状態次元 `r` と推論ツール**: `r = 2` 固定（論文にデータ駆動の選択規則なし）。本パッケージで唯一、専用推論関数・`conformal_inference()` とも未対応
- **Algorithm 7 厳密 M-step と全観測セル EM**:
  - E-step は Joseph form のカルマンフィルタ + RTS スムーザ。M-step は平滑化事後共分散 `P^s`・`P_cross` を含む二次モーメントで `W`・`mu`・`[A C]`・対角 `Q`/`R` を更新し、RTS 1 歩戻しで `m0`・`P0` も更新する（`P^s` を落とすと `Q` が下限へ潰れる）。欠測セルは観測された期間のみで更新する（Shumway & Stoffer 1982）
  - ドリフト `C` と全観測セル EM（論文基本形 `EM_pre` でなくマスクされた全観測使用）を維持する（staggered の単一フィット性・高水準パネルでの精度のため）
- **ユニット切片 `mu_i`（`intercept = TRUE` 既定）**: 切片なしだと水準を `z0`・`C`・単位根近傍 `A` で運び合い尤度が平坦化して停止時点依存になるため、`[m_t, 1]` 回帰で `mu_i` を同時推定する（`intercept = FALSE` で論文基本形）
- **収束判定**: `em_tol`（既定 `1e-6`）は当てはめ値 `mu + W z` の反復間変化を**ユニット平均まわりの変動ノルム**で割って判定する（生ノルムで割ると高水準パネルで早期停止する）。`em_iter`（既定 `2000`）到達時は警告する
- **カルマンゲイン最適化**: `kalman_smoother_cpp` は対角 `R` に対して push-through 恒等式で `r × r` 逆行列に落とし、EM ループ内では `R` 対角成分をベクトルのまま持ち回る（返却時のみ行列化、E-step には `Y_obs - mu` を渡す）
- **`Y_hat` 規約**: TASC のみ `Y_hat` が `T × N` の**全ユニット**行列（`mu` 込み）を返すため、消費側（`synthetic_outcomes()`・`export_json()`）が `idx_tr` 列を選択する
- **staggered と初期化**: 観測マスク `Y_obs`（採用時点以降 `NA`）による単一フィットを `.impute_cohort_fits()` で分割し Clarke 加重で集計する（`control_group` 引数なし）。SVD / `mu` 初期化はグローバル `T_pre`（最早採用時点 - 1）窓を使い、never-treated が `r` 以上あればそれのみ、不足すれば全ユニットを使う（未処置期間ゼロのユニットのみ `.drop_unidentified()` で除外し、状態方程式が外挿するため期間は落とさない）
