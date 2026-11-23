---
trigger: glob
globs: "R/si.R,src/si.cpp"
---

# SI（Agarwal et al. 2025）

$$Y_{ti}^{(d)} = \sum_\ell u_{t\ell} v_{i\ell} \lambda_{d\ell} + \varepsilon_{ti}^{(d)}$$

- **既知の制限**: 共変量は sharp・staggered とも未対応（将来実装する場合は原著 Appendix E の時不変共変量を事前アウトカムの下に積む形式 `Z_i = [Y_pre,i; X_i]` になる）
- **Multi-arm (`max(d) > 1`) と staggered 設計**: ユニット因子 $v_{i\ell}$ の arm 非依存性を使い、コントロール arm（`d = 0`）の SVD 基底を全処置 arm で共有する。staggered では `Y_pre_co_g`（`T_pre_g × N_co_g`）のみで SVD を取り（future-treated も `T_pre_g` 以前は未処置のため汚染なし）、`(g, d)` セルの二重加重平均で集計する:
  $$\hat\tau = \frac{\sum_g\sum_d N_{tr,g,d} T_{post,g} \hat\tau_{g,d}}{\sum_g\sum_d N_{tr,g,d} T_{post,g}}$$
- **特異値ランク `k` の選択（Gavish–Donoho 2014 SVHT）**: SVD 1 回で済む `tau = omega(beta) * median(singular values)`（Eq. 5 近似）を採用。staggered ではコホートごとの `Y_pre_co_g` に適用し `min(T_pre_g, N_co_g) - 1` で上限を切るため、**staggered には単一の `k` が存在せず `fit$k` はコホート名付き整数ベクトル**（`.si_k_by_cohort()`）になる
- **`si_inference()`（`.unit_inference()` 共通ドライバ）**:
  - 処置（multi-arm は arm ごとに層別抽出）と対照を再抽出し、`N_tr = 1`（いずれかの arm が 1 でも）はエラー（sharp なら `conformal_inference()` へ誘導）、staggered jackknife は全コホート `N_tr,g >= 2` 必須
  - SVD 基底の shared donor pool 前提を保つため、bootstrap はコントロール列を全 arm・全コホートで**共有サンプリング**する。`.si_unit_estimate()` はコホート内の全 arm の処置列を 1 回の `si_pcr_cpp()` でまとめて解く（fit のコホート別 `k` を固定使用）
