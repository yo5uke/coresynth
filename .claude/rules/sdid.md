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

## bootstrap・jackknife（Arkhangelsky et al. 2021 Algorithms 2–3、Clarke et al. 2024 同）

**1.1.0.9000 で論文どおりに作り直した**。旧実装は bootstrap が対照だけを引いて処置群平均と ζ² を固定し、staggered ではコホートごとに独立に対照を引き直していた。jackknife も対照だけを落とし、しかも ω・λ を推定し直していた。被覆率（名目 95%、2 因子 DGP）は bootstrap が N_tr=5 で 0.74、20 で 0.81、staggered（3 コホート×5）で 0.68、jackknife が 0.83 / 0.90 / 0.79（global 0.89）。現行は 0.93 / 0.97 / 0.94、jackknife 0.91 / 0.93 / 0.88（global 0.97）（`sandbox/31_inference_coverage_sdid.R`）。GSC・SI と同じ共通ドライバ `.unit_inference()`（`R/unit_inference.R`、`object-model.md`）に載せてある

- **`method = "bootstrap"`**（Algorithm 2）: 全 N ユニットを**処置・対照をまとめて**復元抽出し（`ctx$pooled`。GSC・SI の層別抽出と違うのは論文の字面どおりだから）、処置 0 または対照 2 未満の標本は**引き直す**（論文の "go to 2"）。`.sdid_unit_estimate()` が共変量の partial-out（fit が共変量を持つとき、抽出後のパネルで `.sdid_partial_out_staggered()` から）・σ̂²・ζ²（`scm_fit()` で指定されていなければ）・ω・λ をすべて推定し直す。重複して引かれたユニットは別の列。staggered は引いた処置から採用時点ベクトルとコホートを組み直し（Clarke Algorithm 2）、引かれなかったコホートは落ちる
- **`method = "jackknife"`**（Algorithm 3）: 処置も対照も 1 つずつ落とし、**ω（残ったドナーで再正規化）と λ は fit の値で固定**して重み付き DiD だけ計算し直す（`.sdid_jack_estimate()`、`sdid_estimate_cpp`）。共変量の β も固定。staggered の `"jackknife"` はコホート別＋デルタ法、`"jackknife_global"` が Clarke Algorithm 3 そのもの（全ユニットを全コホートから同時に落とす）
- **N_tr=1 は bootstrap・jackknife とも明示エラー**で `"placebo"` へ誘導（論文: jackknife は N_tr=1 で未定義、bootstrap は多数の処置ユニットが前提。Stata `sdid` も bootstrap は処置 2 以上を要求）。**staggered の jackknife はどのコホートにも処置 2 以上が必要**（Clarke Algorithm 3。Stata も同じ制約）— この規則は GSC・SI の jackknife にも当てている
- fit は再推定に使う `zeta2_user`・`zeta_t` を持つ（staggered も）。共変量の経路でしか使わない `Y_raw`・`X_arr` は、共変量なしの fit では `NULL`
- **`sdid_event_study()`**（sharp 限定）: 同じ bootstrap をユニットまたは**クラスタ単位**で行い、期間別効果も返す。こちらは処置 0・対照 2 未満の draw を引き直さずに捨てて `n_boot_dropped` に数える（`sdid_inference()` は引き直す）。各 draw は「抽出したデータフレームに `scm_fit()` をかけ直した結果」と丸め誤差の範囲で一致する（テストで照合済み）。期間別効果は λ 整列済みのギャップで、事後平均が τ̂ に一致する。`sdid_inference()` は `cluster` を受け付けず、渡すと `sdid_event_study()` へ誘導するエラー
- クラスタの抽出順序は `levels(factor(列))`。`split(data, data[[cluster]])` と同じ順序なので、同じ seed なら手書きのループと同じ draw 列になる。数値の ID を文字列のまま `sort()` すると `"10" < "2"` の順になり、ずれる

**placebo は今回の対象外**（Algorithm 4 との差は backlog）: λ・ζ² を fit の値で固定し、N_tr>1 でも疑似処置を 1 ユニットにしている。実測で N_tr=5・20 の被覆が 0.97〜1.00（SE が SD の 1.3〜2.6 倍）と保守的、N_tr=1 では 0.81（SE 0.33 vs SD 0.49）と逆に過小

## staggered placebo

never-treated 集合は全コホートの `idx_co` の共通部分（既定の `"never_treated"` では全コホート共通。非推奨の `"clean"` では future-treated が後期コホートの `idx_co` に含まれないため自然に除外される）。各コホートで `sdid_placebo_cpp`（λ・ζ² 固定）を呼び、加重平均でグローバル placebo ATT を作る。p 値はこの置換分布から計算する。OpenMP 並列化の設計は `scm.md` の placebo 節と同じパターンを共有している。

## plot 固有の実装

- **weights の 2 パネル表示**: `coord_flip` は `facet_wrap(scales="free")` と併用できないため、SDID の weights プロットだけ native 横棒で ω/λ を別々に描画している
- **`align`**: SDID の推定量は定義上「整列した状態での post 平均差」（Arkhangelsky et al. 2021）なので、λ 加重で整列すると post 平均ギャップが τ̂ と厳密に一致する。他手法の単純プリ平均差とは意味が異なる。一般的な plot 規約は `plot.md` 参照
