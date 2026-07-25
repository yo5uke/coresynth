# 理論的整合性（論文照合済み）

全実装ファイル（R + C++）と対応論文の照合結果。全手法横断の一覧（各手法の理論的注意点・個別の設計判断は `.claude/rules/<method>.md` 参照）。

| 手法 | 検証項目 | 根拠 |
|---|---|---|
| SDID | σ̂²（Eq. 2.2）・ζ²・QP 正則化・推定量（Eq. 2.4–2.5） | Arkhangelsky et al. (2021) |
| SCM | 座標降下 11 点グリッド・OOS V 選択・低ランク FISTA パス・`proj_simplex` 単調早期 break | Abadie et al. (2010) |
| GSC | EM ループ（E: truncated SVD、M: ridge OLS）・処置ローディング OLS・反事実構成 | Xu (2017) |
| SI-PCR | `Y_pre_co^T` の SVD から疑似逆行列構成（`Pinv = U_k diag(1/S_k) V_k^T`） | Agarwal et al. (2025) §3.3 |
| MC | Soft-Impute（fill-in + SVD 閾値）アルゴリズム | Mazumder et al. (2010) |
| TASC（C++） | `kalman_smoother_cpp`（Joseph form）・`arma::find_finite` による NA 除外 | Shumway & Stoffer (2000) |
| SCM-Design | base/weakly_targeted/unit_level 3 バリアント（eq.7/9/10 を S 固定下で厳密最適化）・順列検定・共形 CI | Abadie & Zhao (2026) + 付属コード SCDesign |
| Staggered 加重 | `Σ N_tr,g·T_post,g·τ̂_g / Σ N_tr,g·T_post,g` を全手法で正しく実装 | Clarke et al. (2023) |
| SDID 推論（sharp） | placebo（時間重み固定・ユニット重み再計算）・bootstrap・jackknife | Clarke et al. (2023) |
| SDID 推論（staggered） | bootstrap/jackknife/jackknife_global + placebo（never-treated intersection + per-cohort 置換検定） | Clarke et al. (2023) |
| SI-PCR 多アーム（sharp） | ユニット因子 v_il の arm 非依存性: コントロール arm SVD 基底を全 arm で共有。加重平均 ATT | Agarwal et al. (2025) §5 |
| SI-PCR 多アーム + staggered | (g,d) セルの二重加重平均・コホート内 SVD 基底共有・shared donor pool を保つ推論 | Agarwal et al. (2025) §5 + Clarke et al. (2023) |
| Conformal | 全 T 再推定・moving-block 置換・test inversion CI | Chernozhukov, Wüthrich & Zhu (2021) |
| SCM 検証ツール | in-time placebo（backdating）・ドナー LOO（V 固定, fn.20）・placebo gap プロット（Fig 4–8, mspe_prune 20/5/2） | ADH (2015)・ADH (2010) §3.4 |
| SCM staggered | partially pooled SCM（Eq. 6 の ν 目的関数）・intercept shift（Eq. 8–10）・wild bootstrap（Eq. 12–13）・コホート内完全プーリングの正当性（App A.2）。原著者実装 augsynth と重み・ν̂ が一致 | Ben-Michael, Feller & Rothstein (2022) + 付属コード augsynth |
