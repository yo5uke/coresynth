---
paths:
  - "R/gsc.R"
  - "src/gsc.cpp"
---

# GSC（Xu 2017）

$$Y_{it} = \delta_{it} D_{it} + x_{it}'\beta + \lambda_i' f_t + \varepsilon_{it}$$

共変量なしは truncated SVD 3 ステップ、共変量ありは EM ループ（E-step: SVD、M-step: リッジ OLS）。**因子数 r は既定で CV 選択**（`r = NULL`）。推論は `gsc_boot()`（sharp のみ、H0 下の parametric bootstrap）と `gsc_inference(method="bootstrap"/"jackknife"/"jackknife_global")`（sharp + staggered、非パラメトリック）。

**Staggered 設計**: `Y_co_g` は全 T 期間を使う（future-treated の処置後汚染を SDID と同様に許容する設計）。共変量は全ユニット T×N×p 配列を per-cohort にスライスして EM に渡す。

**因子正規化の注意**: 実装は `F = U_k · diag(D_k)` で Xu (2017) の識別制約（`F'F = T·I_r`）を満たさないが、反事実推定は F の列空間への射影であり正規化定数に依存しないため問題ない。

**因子数の CV**（`.gsc_cv_r()`、Xu 2017 Algorithm 1）: 候補 r ごとに ①コントロールで IFE を推定 ②プリ期間 s を 1 つずつ hold out し残りで処置ユニットの loading を OLS 再推定 ③s を予測して誤差を貯める、で MSPE(r) = Σ_s Σ_i e²_is / T_pre を最小化する r を選ぶ（分母が T_pre のみなのは原論文 Step 3 通り。N_tr は候補間で一定なので argmin に影響しない）。**候補は 1 以上** — r=0（因子なし＝2-way FE）は `gsc_ife_cpp` が rank-r 切断 SVD を取る構造上表現できない。上限は `min(5, T_pre-1, N_co)`（1 期間 hold out すると loading の OLS に使えるのは T_pre-1 期間）。staggered でも**グローバルなプリ期間で 1 回だけ**回す — 全処置ユニットがその窓では未処置なので定義が一意に定まる。選択表は `fit$r_cv`。既知ランクの DGP で r=1〜4 を完全復元することを確認済み
