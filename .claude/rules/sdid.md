---
paths:
  - "R/sdid.R"
  - "src/sdid.cpp"
---

# SDID（Arkhangelsky et al. 2021）

$$\hat\tau^{sdid} = (\bar{Y}_{tr,post} - \hat\lambda^\top Y_{tr,pre}) - (\hat\omega^\top \bar{Y}_{co,post} - \hat\omega^\top Y_{co,pre}^\top \hat\lambda)$$

ユニット重み（Eq. 2.1）は列方向デミーニング QP。推論は `sdid_inference(method="placebo"/"bootstrap"/"jackknife"/"jackknife_global")` で sharp/staggered とも全 4 手法に対応。**Covariates**: Clarke et al. (2023) の global `Wit=0` partial-out（sharp = コントロール×事前期間 OLS、staggered = 全 `Wit=0` 観測での global OLS。SCM staggered も同じ実装を再利用している）。partial-out には数値安定化用の微小 ridge（`diag(1e-8)`）を加えている（論文は plain OLS、推定への影響は誤差範囲）。

## staggered placebo

never-treated 集合は全コホートの `idx_co` の共通部分（`control_group="clean"` では future-treated が後期コホートの `idx_co` に含まれないため自然に除外される）。各コホートで `sdid_placebo_cpp`（λ・ζ² 固定）を呼び、加重平均でグローバル placebo ATT を作る。p 値はこの置換分布から計算する。OpenMP 並列化の設計は `scm.md` の placebo 節と同じパターンを共有している。

## plot 固有の実装

- **weights の 2 パネル表示**: `coord_flip` は `facet_wrap(scales="free")` と併用できないため、SDID の weights プロットだけ native 横棒で ω/λ を別々に描画している
- **`align`**: SDID の推定量は定義上「整列した状態での post 平均差」（Arkhangelsky et al. 2021）なので、λ 加重で整列すると post 平均ギャップが τ̂ と厳密に一致する。他手法の単純プリ平均差とは意味が異なる。一般的な plot 規約は `plot.md` 参照
