---
paths:
  - "R/si.R"
  - "src/si.cpp"
---

# SI（Agarwal et al. 2025）

$$Y_{ti}^{(d)} = \sum_\ell u_{t\ell} v_{i\ell} \lambda_{d\ell} + \varepsilon_{ti}^{(d)}$$

SI-PCR: $Y_{pre,co}^\top$ の truncated SVD → 疑似逆行列でドナー重みを計算する。デフォルト $k = \lfloor\sqrt{\min(T_{pre}, N_{co})}\rfloor$。推論は `si_inference(method="bootstrap"/"jackknife"/"jackknife_global")` — sharp・staggered・multi-arm・staggered×multi-arm すべてに対応。

**既知の制限**: 共変量は未対応（sharp/staggered とも。Agarwal et al. 2025 の原著が covariate なしの定式化のため）。

**Multi-arm (K>1)**: ユニット因子 $v_{i\ell}$ が arm に依存しないという性質を使い、コントロール arm (d=0) の SVD 基底を全 K 処置 arm で共有する。`d` が多値（`max(d) > 1`）なら自動的に multi-arm 経路へ dispatch する。

**Staggered 設計**: `Y_pre_co_g`（T_pre_g × N_co_g）のみで SVD を取る — future-treated も T_pre_g 以前は未処置なので汚染がない。**Staggered × multi-arm** はコホート g × アーム d の (g,d) セルごとに SI-PCR を行い、コホート内で SVD 基底を共有する。集計式:

$$\hat\tau = \frac{\sum_g\sum_d N_{tr,g,d} T_{post,g} \hat\tau_{g,d}}{\sum_g\sum_d N_{tr,g,d} T_{post,g}}$$

**Multi-arm 推論の注意**: bootstrap はコントロール列を全 arm（staggered ではコホート内全 arm）で**共有サンプリング**する必要がある — SVD 基底が shared donor pool を前提にしているため、arm ごとに独立リサンプルすると基底の整合性が崩れる。sharp multi-arm には `jackknife_global` が無い（staggered 専用のため）。

**k の選択**: Agarwal et al. は CV・universal SVHT・elbow の 3 つをデータ駆動の方法として挙げる。採用したのは 2 番目 — `tau = omega(beta) * median(singular values)`、`beta = min(dim)/max(dim)`、`omega` は原論文 Eq.5 の閉形式近似（beta ∈ (0,1] で誤差 0.0001 以下と報告されている）。**SVD 1 回で済む**のが決め手（CV だと候補ごとに PCR フィットが要る）。旧既定 `floor(sqrt(min(T_pre, N_co)))` は上記いずれでもない独自式だった。既知ランクの DGP で r=1〜4 を完全復元（旧式は常に 4 を返していた）。staggered/multi-arm はコホートの `Y_pre_co_g` に対して同じ閾値を当て、`min(T_pre_g, N_co_g) - 1` で上限を切る
