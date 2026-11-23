---
trigger: glob
globs: "R/sdid.R,src/sdid.cpp,src/optim.cpp"
---

# SDID（Arkhangelsky et al. 2021; Clarke et al. 2024）

$$\hat\tau^{sdid} = (\bar{Y}_{tr,post} - \hat\lambda^\top Y_{tr,pre}) - (\hat\omega^\top \bar{Y}_{co,post} - \hat\omega^\top Y_{co,pre}^\top \hat\lambda)$$

## 重み推定と数値安定性

- **二方向 FE 付き partial-out（Kranz 2022 / Clarke et al. 2024 S.2.2）**: `covariates` は sharp = コントロール×事前、staggered = 全 `Wit = 0` 観測上で `.twfe_partial_out()`（非平衡セルのため交互射影、`diag(1e-8)` ridge 付き、FE に吸収される共変量は `beta = 0` にして警告）により除去する。SCM staggered も共用
- **行平均デミーニングと横長最適化（`solve_simplex_ridge_ls()` / `sdid_placebo_cpp()`）**: `1'w = 1` より `(A - a1')w - (b - a) = Aw - b` は厳密に不変。行平均を引かずに解くとユニット水準や共通トレンドが Lipschitz 定数を膨らませ FISTA が tol `1e-6` で早期停止するため必須。また横長 `A`（ユニット重み）は N×N Gram を作らず勾配 `A'(Ay) + rho*y - c` で回す
- **対照 1 つ以上で動作**: ドナー 1 つでも `omega = 1` の `lambda` 加重 DiD として定義される（SCM/GSC/SI の最小対照数 2 とは異なる）。`unit_weights` は他手法と同じ名前つき数値ベクトル（`as.numeric()` 正規化）で返す
- **`synthdid` パッケージとの意図的差異**: 推定量の定義は同一だが、coresynth は FISTA で厳密解まで解く（synthdid の Frank-Wolfe + sparsify は近似誤差が残る）、`sigma_hat^2` 分母は `n`（synthdid は `n - 1`）、推論では `sigma_hat^2` / `zeta^2` も再推定する（`zeta2` 指定時除く）、`N_tr = 1` の placebo は全対照列挙で乱数不使用

## 推論契約（`.unit_inference()`・`sdid_event_study()`）

- **`method = "bootstrap"`（Algorithm 2）**: `ctx$pooled` で全 N ユニットをまとめて復元抽出し、処置 0 または対照 0 の標本は引き直す。各 draw で partial-out（共変量あり時）・`sigma_hat^2`・`zeta^2`・`omega`・`lambda` をすべて再推定し、staggered は引いた処置から採用時点とコホートを組み直す
- **`method = "jackknife"`（Algorithm 3）**: `omega`（残余ドナーで再正規化）・`lambda`・共変量 `beta` を fit 値に固定して重み付き DiD のみ再計算する（`.sdid_jack_estimate()`）。`N_tr = 1` は placebo へ誘導してエラー、staggered は全コホート `N_tr,g >= 2` 必須（`"jackknife"` はコホート別 + デルタ法、`"jackknife_global"` は Clarke Algorithm 3）
- **`method = "placebo"`（Algorithm 4）**: 対照のみから `N_tr` 個を疑似処置に割り当て推定量全体を再推定する（staggered は処置ユニットの採用時点をランダムに引き継ぐ、`N_co > N_tr` 必須）。`N_tr = 1` は全対照列挙で名前つき、`N_tr > 1` は `n_boot` 組の非復元抽出。p 値は `(1 + #{|tau_b| >= |tau_hat|}) / (B + 1)`、SE は `1/B` 正規化。旧 `sdid_placebo_cpp()` は export 互換のため残置しているが `sdid_inference()` では使わない
- **共変量と補助データ**: fit は再推定用に `zeta2_user`・`zeta_t` を保持し、`Y_raw`・`X_arr` は共変量ありのときのみ保持する（なしでは `NULL`）
- **`sdid_event_study()`（sharp 限定）**: 処置 0 / 対照 0 の draw は引き直さず `n_boot_dropped` に数える。クラスタ順は `split()` と一致させるため `levels(factor(cluster))` を使う（`sort(as.character())` 禁止）。期間別効果は `lambda` 整列済みギャップで、事後平均が `tau_hat` と一致する
- **plot 固有**: weights の 2 パネル表示は `coord_flip()` と `facet_wrap(scales = "free")` の非互換回避のため native 横棒描画、trend/gap の `align` は `lambda` 加重整列を使う
