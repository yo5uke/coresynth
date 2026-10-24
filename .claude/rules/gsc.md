---
paths:
  - "R/gsc.R"
  - "src/gsc.cpp"
---

# GSC（Xu 2017）

$$Y_{it} = \delta_{it} D_{it} + x_{it}'\beta + \lambda_i' f_t + \varepsilon_{it}$$

共変量なしは truncated SVD 3 ステップ、共変量ありは Bai (2009) の反復（因子 = `Y − Xβ` の切断 SVD、β = `(Σ X'M_F X)⁻¹ Σ X'M_F Y`）。**因子数 r は既定で CV 選択**（`r = NULL`）。推論は `gsc_boot()`（sharp のみ、Xu 2017 Algorithm 2 の parametric bootstrap）と `gsc_inference(method="bootstrap"/"jackknife"/"jackknife_global")`（sharp + staggered、ユニット再抽出）。

**共変量ループの β 更新は集中化形**（1.1.0 で変更）: 旧実装は `β = (X'X)⁻¹X'(Y − FΛ')` の素朴な交互最小二乗で、固定点は同じだが共変量が因子・ローディングと相関すると極端に遅い。`max_iter = 50` で黙って打ち切っていたため、Xu (2017) 式(3) の DGP（論文自身の Monte Carlo）で β̂=(0.789, 2.768)（真値 1, 3）、ATT bias +0.062 だった。集中化形は同じ固定点に 39〜49 反復で着く（ALS は 221 反復）。既定は `max_iter = 1000`・`tol = 1e-8`、戻り値に `n_iter`/`converged`、R 側は `.gsc_warn_unconverged()` で警告する。収束後は fect `method="gsynth", force="none"`（同じ r）と 30 draws で ATT 平均差 0.029（fect の tol=1e-3 の分）

**Step 1 と Step 2 を分けてある**: `.gsc_fit_controls()` が対照だけの IFE（`gsc_ife_cpp` に 0 列の処置を渡す）、`.gsc_counterfactual()` が処置側（C++ の `gsc_treated_cpp`、非 export。`gsc_ife_cpp` 内部の処置ステップと同じ static 関数）。staggered の `"never_treated"` では全コホートの対照が同一なので、IFE を 1 回だけ当ててコホートごとにローディングだけ求める（コホート数倍の削減、結果は不変）。`"clean"` はコホートごとに対照が違うので従来どおり個別に当てる

**Staggered 設計**: 因子は対照の全 T 期間から、処置ローディングは各コホートの事前窓から。共変量は全ユニット T×N×p 配列をスライスして渡す。fit は `control_group` を記録する（推論の再推定が同じ規則を使うため）

**`gsc_boot()` は Algorithm 2 どおり**（1.1.0 で作り直し）: 旧実装は iid ガウス・共通 σ̂ で、系列相関 ρ=0.5 のとき被覆率 77%（N_tr=1）。論文は ①対照を 1 つ抜いて「処置扱い」にし、残りの対照から N_co 個を復元抽出して当てた予測誤差ベクトル（B 本）を処置側に、②対照の in-sample 残差ベクトルを対照側に、**時系列ごと**再抽出する（系列相関を保つ、処置側は予測誤差なので分散が大きい）。試作で被覆 91%/94%（ρ=0/0.5）。CI は `ATT_obs − q(1−α/2)`〜`ATT_obs − q(α/2)`。論文の字面は ATT＋分位点の percentile 法だが、複製は効果ゼロで生成するので ATT* は推定誤差の分布であり、basic 形の方が枢軸的に正しい（対称なら一致）。分位点そのもの（0 中心の受容域）を CI として返していた不具合を 1.0.0 で修正済み — 戻さないこと

**`gsc_inference()` は処置も再抽出する**（1.1.0 で変更、共通ドライバ `.unit_inference()` は `R/unit_inference.R`）: 旧実装は対照だけを引き処置を固定していたため、処置側のノイズ（事後ショック・事前期間から推定するローディングの誤差）が SE に入らず、N_tr=1 で被覆 30%、staggered（3 コホート×5）で 42%。bootstrap は処置と対照を層別に復元抽出し、staggered では処置を全コホートからプール（Clarke Algorithm 2・fect と同じ。コホートは引いた処置から組み直し、重みも再計算）、対照は 1 回だけ引いて全コホートで共有する（旧実装はコホートごとに独立に引き、共有ドナー由来のコホート間相関を落としていた）。jackknife は対照も処置も 1 つずつ落とす。**N_tr=1 は明示エラー**で `gsc_boot()`/`conformal_inference()` へ誘導する（Xu S.3.2: ユニット再抽出が妥当なのは N_tr が大きいとき）。staggered の jackknife は、どのコホートにも処置 2 以上が必要（Clarke et al. 2024 Algorithm 3、処置 1 つのコホートは落とすと効果が未定義になる。SDID と共通）。層別 bootstrap の被覆は N_tr=5 で 0.88、20 で 0.95

**`conformal_inference()` は共変量込みで再推定する**（1.1.0 で修正）: 旧実装は β を落として IFE を当て直しており、共変量つき fit で再推定の事前 RMSE 2.716（fit は 0.468）、CI は共変量なしと同一になっていた。処置平均の系列には共変量の処置平均を当てる

**加法的固定効果 `force`**（1.1.0 で追加、Xu 2017 Remark 3 の拡張モデル `Y = δD + x'β + μ + α_i + ξ_t + λ'f + ε`）: `"none"`（既定・後方互換）/ `"unit"` / `"time"` / `"two-way"`。対照パネルを within 変換（`"unit"` は各ユニットの時間平均、`"time"` は各期の横断平均、`"two-way"` は両方を引く）してから同じ Bai 反復を回す — 平衡パネルなので FWL により FE＋因子の同時最小二乗と厳密に一致し、因子は自動的に定数と直交する。μ・α_i・ξ_t は収束後の残差 `Y − Xβ − FΛ'` の平均から戻す。処置ユニットは事前期間で `Y − Xβ − μ − ξ_t` を `[1, F̃]`（ユニット効果がなければ `F̃` のみ）に回帰して (α_i, λ_i) を得る。以前の R 側 within 変換の試作が gsynth を再現できなかったのは、処置ユニットを全期間でデミーンしていた（事後期間が混ざる）ため。`force ≠ "none"` では r=0 も許し、CV 候補も 0 から（gsynth の `r = c(0, 5)` と同じ）。**gsynth/fect と r=0,1,2・`"unit"`/`"time"`/`"two-way"` の全組み合わせで ATT が 1e-15 まで一致**（2 因子＋ユニット・時間効果の DGP）。**`"none"` だけは gsynth と違う**: gsynth の `"none"` は総平均 μ を推定しており、Y から対照の総平均を引いてから coresynth `"none"` を当てると機械精度で一致する。coresynth の `"none"` は Xu 式(1) どおり μ も持たない。fit は `force`・`mu`・`alpha_co`・`alpha_tr`・`xi` を持ち、`gsc_boot()`・`gsc_inference()`・`conformal_inference()`・staggered の各再推定は同じ `force` を使う

**因子構造が弱いパネルでは r と force に強く依存する**: Prop 99 では gsynth 自身が force 次第で +0.549（unit）/ −1.873（two-way）/ −20.889（none, r=2）と振れ、coresynth（`"none"`）は CV r=4 で −7.152。この乖離は実装差ではない（同じ r で前処置 RMSPE は coresynth が全 r で gsynth の `"none"` を下回る。r=4 で 1.38 vs 2.30。gsynth の `"none"` は μ 込みなので厳密には別モデル）。GSC が Prop 99 に向かないだけ

**因子正規化の注意**: 実装は `F = U_k · diag(D_k)` で Xu (2017) の識別制約（`F'F = T·I_r`）を満たさないが、反事実推定は F の列空間への射影であり正規化定数に依存しないため問題ない。

**因子数の CV**（`.gsc_cv_r()`、Xu 2017 Algorithm 1）: 候補 r ごとに ①コントロールで IFE を推定 ②プリ期間 s を 1 つずつ hold out し残りで処置ユニットの loading を OLS 再推定 ③s を予測して誤差を貯める、で MSPE(r) = Σ_s Σ_i e²_is / T_pre を最小化する r を選ぶ（分母が T_pre のみなのは原論文 Step 3 通り。N_tr は候補間で一定なので argmin に影響しない）。**候補は 1 以上** — r=0（因子なし＝2-way FE）は `gsc_ife_cpp` が rank-r 切断 SVD を取る構造上表現できない。上限は `min(5, T_pre-1, N_co)`（1 期間 hold out すると loading の OLS に使えるのは T_pre-1 期間）。staggered でも**グローバルなプリ期間で 1 回だけ**回す — 全処置ユニットがその窓では未処置なので定義が一意に定まる。選択表は `fit$r_cv`。既知ランクの DGP で r=1〜4 を完全復元することを確認済み
