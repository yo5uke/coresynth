---
paths:
  - "R/sdid.R"
  - "src/sdid.cpp"
---

# SDID（Arkhangelsky et al. 2021）

$$\hat\tau^{sdid} = (\bar{Y}_{tr,post} - \hat\lambda^\top Y_{tr,pre}) - (\hat\omega^\top \bar{Y}_{co,post} - \hat\omega^\top Y_{co,pre}^\top \hat\lambda)$$

ユニット重み（Eq. 2.1）は列方向デミーニング QP。推論は `sdid_inference(method="placebo"/"bootstrap"/"jackknife"/"jackknife_global")` で sharp/staggered とも全 4 手法に対応。**Covariates**: Kranz (2022) / Clarke et al. (2024, S.2.2 "projected") の partial-out（sharp = コントロール×事前期間、staggered = 全 `Wit=0` 観測）。**回帰は `Y = Xβ + μ_i + λ_t` の二方向 FE つき**（`.twfe_partial_out()`、非平衡セルなので交互射影で within 変換）。v1.0.0 以前は FE なしの OLS で、β が共変量と結果の水準・トレンドの相関を吸って β̂=3.75（真値 1）、ATT 2 → SCM 1.65 / SDID 1.18 と崩れていた。FE に吸収される共変量（時間不変・共通時系列）は β=0 にして警告。SCM staggered も同じ実装を再利用している。数値安定化用の微小 ridge（`diag(1e-8)`）あり。

## 重み QP のソルバ

`sdid_unit_weights_cpp`/`sdid_time_weights_cpp` は `solve_simplex_ridge_ls()`（`optim.cpp`）で `min ½‖Aw‖² + ½ρ‖w‖² − (A'b)'w` を解く。A が横長（ユニット重み: T_pre 行 × N_co 列）なら N×N の Gram を作らず勾配 `A'(Ay) + ρy − c` と `λ_max(AA') + ρ` で FISTA を回す — 反復列は密形式と丸め誤差を除いて同一。N_co=300 で placebo/jackknife が 4〜8 倍速。

## bootstrap の 2 系統

- **`sdid_inference(method = "bootstrap")`**: 統制ユニットだけを復元抽出し、処置群平均系列と ζ² は元の推定値に固定する。`cluster` は受け付けず、渡すと `sdid_event_study()` へ誘導するエラーになる
- **`sdid_event_study()`**（sharp 限定）: 処置群と統制群の両方をユニット単位またはクラスタ単位で復元抽出する。重複して引かれたユニットは別の列（新しい ID）として扱い、σ̂²・ζ²（`scm_fit()` で指定されていなければ）・ω・λ を毎回推定し直す。処置群が空、または統制が 2 未満の draw は捨てて `n_boot_dropped` に数える。共変量があれば抽出後のパネルで `.sdid_partial_out_staggered()` からやり直す。各 draw は「抽出したデータフレームに `scm_fit()` をかけ直した結果」と丸め誤差の範囲で一致する（テストで照合済み）。期間別効果は λ 整列済みのギャップで、事後平均が τ̂ に一致する
- クラスタの抽出順序は `levels(factor(列))`。`split(data, data[[cluster]])` と同じ順序なので、同じ seed なら手書きのループと同じ draw 列になる。数値の ID を文字列のまま `sort()` すると `"10" < "2"` の順になり、ずれる
- 共変量の経路でしか使わない `Y_raw`・`X_arr` は、共変量なしの fit では `NULL`

## staggered placebo

never-treated 集合は全コホートの `idx_co` の共通部分（既定の `"never_treated"` では全コホート共通。非推奨の `"clean"` では future-treated が後期コホートの `idx_co` に含まれないため自然に除外される）。各コホートで `sdid_placebo_cpp`（λ・ζ² 固定）を呼び、加重平均でグローバル placebo ATT を作る。p 値はこの置換分布から計算する。OpenMP 並列化の設計は `scm.md` の placebo 節と同じパターンを共有している。

## plot 固有の実装

- **weights の 2 パネル表示**: `coord_flip` は `facet_wrap(scales="free")` と併用できないため、SDID の weights プロットだけ native 横棒で ω/λ を別々に描画している
- **`align`**: SDID の推定量は定義上「整列した状態での post 平均差」（Arkhangelsky et al. 2021）なので、λ 加重で整列すると post 平均ギャップが τ̂ と厳密に一致する。他手法の単純プリ平均差とは意味が異なる。一般的な plot 規約は `plot.md` 参照
