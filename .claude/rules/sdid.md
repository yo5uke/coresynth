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

**対照 1 つでも動く**（2026-11-15）: SDID はドナー 1 つでも定義される（ω = 1、λ 加重の DiD）。以前は sharp で対照 2 未満、staggered で never-treated 2 未満をエラーにし、bootstrap も対照 1 の標本を引き直していたが、論文・synthdid とも条件は「1 以上」で、不必要に厳しかった。SCM/GSC/SI は 2 以上のまま（重みや因子を選ぶ余地が要る）。`unit_weights` は他手法と同じ名前つき数値ベクトルで返す（以前は N×1 行列に names 属性が付いた形だった）

## bootstrap・jackknife（Arkhangelsky et al. 2021 Algorithms 2–3、Clarke et al. 2024 同）

**1.1.0.9000 で論文どおりに作り直した**。旧実装は bootstrap が対照だけを引いて処置群平均と ζ² を固定し、staggered ではコホートごとに独立に対照を引き直していた。jackknife も対照だけを落とし、しかも ω・λ を推定し直していた。被覆率（名目 95%、2 因子 DGP）は bootstrap が N_tr=5 で 0.74、20 で 0.81、staggered（3 コホート×5）で 0.68、jackknife が 0.83 / 0.90 / 0.79（global 0.89）。現行は 0.93 / 0.97 / 0.94、jackknife 0.91 / 0.93 / 0.88（global 0.97）（`sandbox/31_inference_coverage_sdid.R`）。GSC・SI と同じ共通ドライバ `.unit_inference()`（`R/unit_inference.R`、`object-model.md`）に載せてある

- **`method = "bootstrap"`**（Algorithm 2）: 全 N ユニットを**処置・対照をまとめて**復元抽出し（`ctx$pooled`。GSC・SI の層別抽出と違うのは論文の字面どおりだから）、処置 0 または対照 0 の標本は**引き直す**（論文の "go to 2"）。`.sdid_unit_estimate()` が共変量の partial-out（fit が共変量を持つとき、抽出後のパネルで `.sdid_partial_out_staggered()` から）・σ̂²・ζ²（`scm_fit()` で指定されていなければ）・ω・λ をすべて推定し直す。重複して引かれたユニットは別の列。staggered は引いた処置から採用時点ベクトルとコホートを組み直し（Clarke Algorithm 2）、引かれなかったコホートは落ちる
- **`method = "jackknife"`**（Algorithm 3）: 処置も対照も 1 つずつ落とし、**ω（残ったドナーで再正規化）と λ は fit の値で固定**して重み付き DiD だけ計算し直す（`.sdid_jack_estimate()`、`sdid_estimate_cpp`）。共変量の β も固定。staggered の `"jackknife"` はコホート別＋デルタ法、`"jackknife_global"` が Clarke Algorithm 3 そのもの（全ユニットを全コホートから同時に落とす）
- **N_tr=1 は bootstrap・jackknife とも明示エラー**で `"placebo"` へ誘導（論文: jackknife は N_tr=1 で未定義、bootstrap は多数の処置ユニットが前提。Stata `sdid` も bootstrap は処置 2 以上を要求）。**staggered の jackknife はどのコホートにも処置 2 以上が必要**（Clarke Algorithm 3。Stata も同じ制約）— この規則は GSC・SI の jackknife にも当てている
- fit は再推定に使う `zeta2_user`・`zeta_t` を持つ（staggered も）。共変量の経路でしか使わない `Y_raw`・`X_arr` は、共変量なしの fit では `NULL`
- **`sdid_event_study()`**（sharp 限定）: 同じ bootstrap をユニットまたは**クラスタ単位**で行い、期間別効果も返す。こちらは処置 0・対照 0 の draw を引き直さずに捨てて `n_boot_dropped` に数える（`sdid_inference()` は引き直す）。各 draw は「抽出したデータフレームに `scm_fit()` をかけ直した結果」と丸め誤差の範囲で一致する（テストで照合済み）。期間別効果は λ 整列済みのギャップで、事後平均が τ̂ に一致する。`sdid_inference()` は `cluster` を受け付けず、渡すと `sdid_event_study()` へ誘導するエラー
- クラスタの抽出順序は `levels(factor(列))`。`split(data, data[[cluster]])` と同じ順序なので、同じ seed なら手書きのループと同じ draw 列になる。数値の ID を文字列のまま `sort()` すると `"10" < "2"` の順になり、ずれる

## placebo（Arkhangelsky et al. 2021 / Clarke et al. 2024、Algorithm 4）

**2026-11-15 に論文どおりに作り直した**。N_tr 個の対照に placebo を割り当て、**対照だけで推定量全体を当て直す**（`.sdid_unit_estimate()`、bootstrap と同じ再推定: 共変量の partial-out・σ̂²・ζ²（未指定なら）・ω・λ、切片つき）。staggered では placebo ユニットが処置ユニットの採用時点ベクトルをランダムな対応で引き継ぐ（Clarke: "based on the actual treatment structure"）。

- **N_tr = 1 は全対照を順に placebo にする**（列挙）。Algorithm 4 の B → ∞ の極限そのもので乱数を使わない。`placebo_effects` は対照名つき。各 placebo は「処置ユニットを除き、その対照を処置扱いにしたデータに `scm_fit()` を掛け直した結果」と**完全一致**する（テストで照合）
- **N_tr > 1 は `n_boot` 組を非復元抽出**（`seed` で再現）。`placebo_effects` は draw ごとの値で名前なし
- never-treated が N_tr + 1 以上必要（論文・synthdid と同じ N_co > N_tr。placebo を除いたドナーが 1 つでも SDID は定義される。以前の `sdid_placebo_cpp()` は N_co = 2 でも動いたので、それを不必要に狭めない）
- p 値は `(1 + #{|τ_b| ≥ |τ̂|}) / (B + 1)`、SE は placebo 分布の 1/B 正規化の標準偏差、CI は正規近似
- 被覆（名目 95%・τ=0・N=40・200 draws）: N_tr=1 で 0.975（SE 0.453 vs SD 0.462）、N_tr=5 で 0.945（0.225 vs 0.232）、staggered 6 ユニットで 0.96（0.270 vs 0.282）
- 速度は placebo 1 つにつきフル再推定 1 回: N_co=299・N_tr=1 で 1.3 秒、N_tr=5・200 draws で 0.6 秒

**旧実装の 3 つのずれ**: ①疑似処置を常に 1 ユニットにしていた（N_tr>1 では平均の分散より大きい量を測るので SE が SD の 1.3〜2.6 倍、被覆 0.97〜1.00）②λ・ζ² を fit の値で固定していた（N_tr=1 で SE 0.33 vs SD 0.49、被覆 0.81）③**ω の QP を列デミーンなしで解いていた** — 本体は ω0 を集中化するので別の推定量になっており、SDID が差し引くはずのユニット固定効果で placebo SE が動いた（水準 SD 0/5/50 で 0.356/0.430/0.474）。②と③が打ち消し合って被覆がそれらしく見えるケースがあった。

`sdid_placebo_cpp()` は export 済みなので残してあるが、`sdid_inference()` はもう使わない（λ・ζ² 固定の LOO placebo。③だけ直して切片つきにした）。

**ユニット水準が大きいと λ・ω の QP の条件が悪くなる**: 理論上は水準に不変だが、時間重みの QP は期間ごとの横断平均しか引かないので、ユニット水準 50〜500 を足すと FISTA の許容誤差 1e-6 の範囲で SE が 4 桁目（相対 1e-3）で動く。

## synthdid（R パッケージ 0.0.9）との仕様差

2026-11-15 に synthdid のソースを読み、同じパネルで実測した（`sandbox/32_sdid_placebo_vs_synthdid.R` は速度のみ）。**推定量の定義（ζ の式・切片の集中化・DiD の形）は同じ**で、差は解法の精度と推論の細部にある。

- **重み QP の解法と精度**: synthdid は Frank-Wolfe（100 反復 → `sparsify_function`（最大値の 1/4 以下を 0 にして再正規化）→ 再び最大 10000 反復、停止は減少量が `min.decrease²` 以下）。coresynth は FISTA で厳密解まで解く。**真の最適解（`solve_simplex_qp` を tol 1e-15・200 万反復で解いた値）からの距離**は、ω の最大誤差が coresynth 7e-7〜1.6e-6（目的関数差 1e-9〜1e-8）、synthdid 既定 1e-4〜1.8e-3（同 8e-6〜9e-2）、sparsify を外し `min.decrease = 1e-10` にしても 1e-7〜3e-4。点推定の差は 12 パネルで平均 1.5e-3・最大 4.8e-3（Prop 99: −15.6046 対 −15.6038、synthdid を厳密に寄せると −15.6063）。**差は synthdid 側の近似誤差で、coresynth を寄せる理由はない**
- **ノイズ水準**: synthdid は全一階差の `sd()`（n−1 で割る）、coresynth は 1/n（比 0.99943。N0=30・T0=30）。ζ をそろえても点推定の差は変わらない（平均 6.6e-4 のまま）ので、差の原因ではない
- **推論での ζ**: synthdid は placebo・bootstrap・jackknife のいずれも元の fit の ζ（`opts`）を固定し、ω・λ は元の重みを初期値に再推定する（jackknife は重みも固定）。coresynth は σ̂・ζ も再推定する（`scm_fit(zeta2 =)` で固定した場合は固定）。bootstrap の SE への影響は 6 パネルの平均で再推定/固定 = 0.992（最大 3.7%）。placebo の SD は列挙で 0.2650 対 0.2649 と一致し、ユニットごとの値の差（最大 0.02）は重み QP の精度による
- **placebo**: synthdid は N_tr = 1 でも毎回ランダムな順列（`placebo_se(est, 200)` は 0.243、列挙の SD 0.265 に対し Monte Carlo 誤差の範囲）、coresynth は N_tr = 1 を全列挙して乱数を使わない。どちらも N_co > N_tr が条件。synthdid は SE のみ、coresynth は置換 p 値と正規近似 CI も返す
- **jackknife**: 式は同じ（`(n−1)/n Σ(θ_i − θ̄)²`）。同じパネルで SE の比は 1.004（重みの差）。N_tr = 1 は synthdid が `NA`、coresynth は明示エラー
- **bootstrap**: 棄却条件はどちらも処置・対照とも 1 以上（論文 Algorithm 2 の "no treated units or no control units"）。coresynth は 2026-11-15 まで対照 2 以上を要求していた
- **既定の推論**: synthdid の `vcov()` は bootstrap（200 回）、`sdid_inference()` は placebo（N_tr = 1 でも動くため）
- **未比較**: 共変量（synthdid は β を重みと同時に最適化、coresynth は Kranz の二方向 FE partial-out）。staggered（synthdid は持たない）

## plot 固有の実装

- **weights の 2 パネル表示**: `coord_flip` は `facet_wrap(scales="free")` と併用できないため、SDID の weights プロットだけ native 横棒で ω/λ を別々に描画している
- **`align`**: SDID の推定量は定義上「整列した状態での post 平均差」（Arkhangelsky et al. 2021）なので、λ 加重で整列すると post 平均ギャップが τ̂ と厳密に一致する。他手法の単純プリ平均差とは意味が異なる。一般的な plot 規約は `plot.md` 参照
