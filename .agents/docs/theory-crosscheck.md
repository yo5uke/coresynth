# 理論的整合性・検証基準（論文照合アーカイブ）

全実装ファイル（R + C++）と対応論文の照合結果、および他パッケージとの基準値。

## 論文照合一覧

| 手法 | 検証項目 | 根拠 |
|---|---|---|
| SDID | σ̂²（Eq. 2.2）・ζ²・QP 正則化・推定量（Eq. 2.4–2.5） | Arkhangelsky et al. (2021) |
| SCM | 座標降下 11 点グリッド・OOS V 選択・低ランク FISTA パス・`proj_simplex` 単調早期 break | Abadie et al. (2010) |
| GSC | Bai (2009) 反復（集中化 β 更新、収束判定つき）・処置ローディング OLS・反事実構成・加法 FE（Remark 3、within 変換＋残差平均）・Algorithm 2 の parametric bootstrap（残差ベクトル再抽出・LOO 予測誤差） | Xu (2017) |
| SI-PCR | `Y_pre_co^T` の SVD から疑似逆行列構成（`Pinv = U_k diag(1/S_k) V_k^T`） | Agarwal et al. (2025) §3.3 |
| MC | Soft-Impute（fill-in + SVD 閾値）アルゴリズム。λ パスの上端は TWFE 残差の最大特異値（MCPanel `mcnnm_lam_range()` と同じ） | Mazumder et al. (2010)・Athey et al. (2021) |
| TASC | `kalman_smoother_cpp`（Joseph form）・`arma::find_finite` による NA 除外・Algorithm 7 の厳密 M-step（P^s 込み、m0/P0 更新、対角 Q・R）＋欠測セル（観測期間のみ）＋ユニット切片 | Rho et al. (2026)・Shumway & Stoffer (1982) |
| SCM-Design | base/weakly_targeted/unit_level 3 バリアント（eq.7/9/10 を S 固定下で厳密最適化）・順列検定・共形 CI | Abadie & Zhao (2026) + 付属コード SCDesign |
| Staggered 加重 | `Σ N_tr,g·T_post,g·τ̂_g / Σ N_tr,g·T_post,g` を全手法で実装 | Clarke et al. (2024) |
| SDID 推論（sharp） | placebo（Algorithm 4: N_tr 個の対照を疑似処置にして推定量全体を再推定）・bootstrap・jackknife | Arkhangelsky et al. (2021)・Clarke et al. (2024) |
| SDID 推論（staggered） | bootstrap/jackknife/jackknife_global + placebo（never-treated から N_tr 個を引き、処置ユニットの採用時点を割り当てて全コホート再推定） | Clarke et al. (2024) |
| SI-PCR 多アーム（sharp） | ユニット因子 v_il の arm 非依存性: コントロール arm SVD 基底を全 arm で共有。加重平均 ATT | Agarwal et al. (2025) §5 |
| SI-PCR 多アーム + staggered | (g,d) セルの二重加重平均・コホート内 SVD 基底共有・shared donor pool を保つ推論 | Agarwal et al. (2025) §5 + Clarke et al. (2024) |
| Conformal | 全 T 再推定・moving-block 置換・test inversion CI（Algorithm 1 の `p > α`） | Chernozhukov, Wüthrich & Zhu (2021) |
| SCM 検証ツール | in-time placebo（backdating）・ドナー LOO（V 固定, fn.20）・placebo gap プロット（Fig 4–8, mspe_prune 20/5/2） | ADH (2015)・ADH (2010) §3.4 |
| SCM staggered | partially pooled SCM（Eq. 6 の ν 目的関数）・intercept shift（Eq. 8–10）・jackknife（App. A.1）・wild bootstrap（Eq. 12 の全ユニット和 ＋ Otsu & Rai 2017 Eq. 8–9）・コホート内完全プーリング（App A.2 の和パラメータ化）。augsynth と不等コホートサイズ含め 4〜5 桁一致 | Ben-Michael, Feller & Rothstein (2022) 本文 + App A.1/A.2 + augsynth・Otsu & Rai (2017) |
| SCM 共変量 | balance で扱う（ADH 2010 S.2.3 の V 加重予測子、BFR S.5.2 Eq. 11 の `cov_method = "balance"`）。sharp の `covariates`（partial-out）は明示エラー。staggered の `partial_out` は Clarke S.2.2 の借用 | Abadie et al. (2010) S.2.3・BFR (2022) S.4/S.5.2 |

## 他パッケージとの照合（既定同士・再検証の基準値）

正準ケースは California Prop 99。`scm_fit(cigsale ~ treated | state + year, method = "scm")` のみで ATT −19.514（ADH 2010 公表値 ≒ −19.6）、ドナー集合も論文と同一。

| 手法 | 相手 | 結果 |
|---|---|---|
| SCM | Synth 1.1.10 / tidysynth | ドナー集合一致・重み 0.02 以内。共通目的関数（前処置 MSPE）で coresynth 3.077 < Synth 3.209（multistart の効果）。ADH 予測子指定では公表重みを max\|diff\| 0.0046 で再現（Synth は 0.0170） |
| 順列検定 | ADH 2010 Fig. 8 | ADH 予測子指定で p=1/39=0.0256 を再現 |
| SDID | synthdid 0.0.9 | Prop 99 で −15.6046 vs −15.6038、重み相関 0.9998。真の最適解（`solve_simplex_qp` tol 1e-15）からの ω 最大誤差は coresynth 7e-7〜1.6e-6 対 synthdid 1e-4〜1.8e-3（12 パネルで点推定差平均 1.5e-3・最大 4.8e-3）。ζ は synthdid が `sd()`（n-1 除算）、coresynth は 1/n 除算（比 0.99943） |
| Conformal | scinference | p=0.096774 完全一致 |
| GSC | gsynth / fect | `force` = unit/time/two-way は r=0,1,2 で ATT 1e-15 一致。`"none"` は gsynth が総平均 μ を持つ分だけ違う（Y から総平均を引けば一致）。共変量つき（Xu 式 3 の DGP、30 draws）で ATT 平均差 0.029 |
| MC | MCPanel / fect | 2.9507 vs 2.9515 / 2.9138（真値 3）。水準別バイアス検証（N_co=20・40 draws）: FE ありで level=0/20 とも −0.018/−0.020（MCPanel 既定 −0.016 と一致）；FE なしは level=20 で +0.230（MCPanel FE 無効 +0.254）。TWFE 残差 σ_max 化後は level 0〜10000 で RMSE 0.3455 と完全不変 |

## 推論の被覆率・検出力の検証基準

- **SDID（`sandbox/31_inference_coverage_sdid.R`、名目 95%、2 因子 DGP）**:
  - bootstrap: N_tr=5 で 0.93、N_tr=20 で 0.97、staggered（3 コホート×5）で 0.94
  - jackknife: N_tr=5 で 0.91、N_tr=20 で 0.93、staggered で 0.88（`jackknife_global` 0.97）
  - placebo（τ=0・N=40・200 draws）: N_tr=1 で 0.975（SE 0.453 vs SD 0.462）、N_tr=5 で 0.945（0.225 vs 0.232）、staggered 6 ユニットで 0.96（0.270 vs 0.282）
- **GSC / SI（`sandbox/30_inference_coverage_gsc_si.R`）**:
  - `gsc_boot()`（N_tr=1）: 被覆 0.91（ρ=0）/ 0.94（ρ=0.5）
  - `gsc_inference()` 層別 bootstrap: N_tr=5 で 0.88、20 で 0.95
  - `si_inference()` 層別 bootstrap: N_tr=5 で 0.88、20 で 0.93、staggered（3 コホート×5）で 0.84
- **SCM**:
  - `scm_rank_pval()`: τ=0 で名目 10%/5% に対しサイズ 7.0%/4.5%（3 ユニット）・10.0%/3.0%（6 ユニット）、検出力 τ=1.5 で 77%/85%
  - `scm_inference(method = "wild_bootstrap")`: 9 設計で 93〜99.5%（J=2 で 93%）
  - `scm_fit(nu = "auto", fixedeff = TRUE)`（BFR S.6 選択あり DGP、N=49・T=39・処置約 20）: two-way FE で separate bias −0.0068/MAD 0.0413 → auto −0.0039/0.0395、factor で ν=1 MAD 0.1492 → auto 0.0607
  - 選択なし DGP での ν̂ 挙動: 独立誤差により ν̂ ≈ 1/√J に収束する（0.733/0.569/0.531 対 1/√J = 0.707/0.535/0.444）。これは既知の性質であり、`scm_balance_frontier()` で確認する
- **TASC**:
  - 水準別 RMSE（20〜30 draws・真値 3、水準 0 / 共通水準 100 / ユニット別水準 100）: `intercept = TRUE` で 0.183 / 0.183 / 0.209。never-treated 0 の staggered パネルでも真値 2 に対し平均 1.93・RMSE 0.27（20 draws）
- **MC**:
  - never-treated 0 の staggered パネル（N=30・T=40・採用 21〜38・真値 2・20 draws）で `.drop_unidentified()` 適用後 ATT 平均 2.00（RMSE 0.20）
