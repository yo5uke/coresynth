---
paths:
  - "R/mc.R"
  - "src/mc.cpp"
---

# MC（Athey et al. 2021, JASA; Mazumder et al. 2010, JMLR）

$$\min_{L, a, b} \frac{1}{2} \|\mathcal{O} \odot (Y - L - a \mathbf{1}_N^\top - \mathbf{1}_T b^\top)\|_F^2 + \lambda_L \|L\|_*$$

反復 SVD ソフト閾値化（Soft-Impute; Mazumder, Hastie & Tibshirani 2010, JMLR）。$\lambda_L$ は $|\mathcal{O}|$ で正規化されない（Athey et al. 2021 式 (4.3) の $\lambda |\mathcal{O}| / 2$ を直接 $\lambda_L$ としてパラメータ化）。**$\lambda$ は既定で K-fold CV 選択**（`lambda = NULL`）、**ユニット・時間の固定効果を罰則外で推定**（`fixedeff = TRUE` 既定）。推論関数は無く、`conformal_inference()` のみ対応。

**観測マスク**: `O` は処置後セルと NA セルの**両方**を 0 にする（`O[is.na(Y)] <- 0`）。NA を観測値 0 として soft-impute に入れると ATT が大きく歪む。`conformal_inference()` の再フィット（`.conformal_refit.coresynth_mc`）でも `M` 内の NA を `O = 0` にマスクする。

**staggered はマスク経由**（TASC と共通）: `O` が各処置ユニットを自身の採用時点から 0 にするため、単一の分解がそのまま staggered fit になる。コホート別の再フィットはしない — `.impute_cohort_fits()`（`utils.R`）が単一 fit の `gap` を採用時点で分割してコホート表を作るだけ。したがって全コホートが同じドナー集合（never-treated ＋ 他処置ユニットの採用前セル）を共有し、`control_group` 引数は持たない（渡すとエラー）。集計は他手法と同じ Clarke 加重。**ユニット等加重にしてはいけない** — コホートの事後長が異なると同一データで他手法と別の推定量になる。

**λ の CV**（`.mc_cv_lambda()`、Athey et al. 2021 S.4）: 観測セル O から K 個の部分集合 O_k を `|O_k| = floor(|O|²/NT)` で引く（fold 内の観測率を元パネルと揃えるため）。各 fold で λ パスを解き O \ O_k で二乗誤差を評価、平均最小の λ を採る。パスは σ_max から 10^-3 倍まで等比 20 点（論文の λ_L=0 は補間解で CV が選ぶことはなく最も重いので置かない）。**論文の warm-start は不採用** — 高速化であって推定量の定義ではなく、`soft_impute_cpp` が常にゼロ初期化のため。**自前の RNG ストリーム**を使い呼び出し側の `.Random.seed` を復元する（λ が呼び出し側の乱数状態に依存すると再現しないため）。fold が成立しない小さすぎるパネルでは旧来の `0.01 × σ_max` にフォールバックし `fit$lambda_cv = FALSE` で告げる

**二元固定効果**（`fixedeff = TRUE` 既定、`soft_impute_fe_cpp`）: JASA 出版版（Athey et al. 2021）では **S.4 式 (4.3) の主推定量そのものが非罰則の二元固定効果を含む**（2018 NBER WP 版では S.4 が `Y = L + e` で FE は S.8.1 の拡張だった。出版版でも S.8.1 に "we regularize the estimates of L*, but do not wish to regularize the estimates of the fixed effects" の議論が残る）。**著者自身の参照実装 MCPanel も `to_estimate_u = to_estimate_v = 1L` が既定**で FE を推定する。実装は FE のバックフィット更新（行・列和の事前計算＋行列ベクトル積 `O * b`, `O' * a`）と soft-impute の SVD ステップを**同一ループ内で交互**に回す（外側で soft-impute を回し直すと CV が現実的な時間で終わらない）。マスクがあるため 2 つの効果は直交せず、1 回の閉形式では解けないので 1 反復あたり 5 スイープ。FE 同士は定数シフト分だけ識別不能だが、完成行列に入るのは和なので影響しない。`fixedeff = FALSE` で FE 無しの純粋低ランクモデル（Mazumder et al. 2010 式 (2)）に戻る。

**旧「ATT 過大推定」の正体は FE 欠落だった**（2026-09-20 に解明）。以前 `mc.md` は「核ノルム正則化の性質、実装の不具合ではない」と記していたが誤り。水準が L に入ると核ノルムがそれを 0 方向へ縮め、水準が正なら ATT が上振れする。N_co=20・40 draws の実測（真値 τ=2）:

| | level=0 | level=20 |
|---|---|---|
| 現行相当（FE 無し） | bias −0.007 | **+0.230 (z=4.5)** |
| MCPanel `to_estimate_u/v=0` | −0.019 | **+0.254 (z=4.8)** |
| FE 有り | −0.018 | **−0.020 (z=−0.4)** |
| MCPanel 既定 | −0.016 | −0.016 |

FE 無し版は MCPanel の FE 無効設定と、FE 有り版は MCPanel 既定と一致する。**以前「CV では直らない」と書いた bias +0.046 は SE≈0.036 で z≈1.3、有意ではなかった**（平均ゼロ DGP だったため水準バイアスが出ていなかった）。

**収束判定は L 単独の相対変化のまま**: 完成行列 M = L + FE に替えると、M はレベルを含んで相対変化が小さく見えるため早く止まり精度が落ちる（tol=1e-12 の解からのずれ 6.9e-5 → 4.7e-4）。λ が大きく L=0 でも FE は exact TWFE と 1e-13 まで一致しており、L 単独で FE が未収束になる問題は無い。

**Choi & Yuan (2024) は不採用**: 脱バイアスは式 (3.1) `M̂ = P_r(Ω^c∘M̃ + Ω∘Y)`（観測セルを戻してランク r 射影）。単独では level=20 で bias +0.123 (z=2.5) と**半減どまりで解消しない**。FE と併用しても FE 単独（−0.020）に対し −0.019 で改善が無い。同論文の本体は MNAR・弱信号下の**推論枠組み**（グループ分割 + 漸近正規性）であり、水準バイアスとは別問題。
