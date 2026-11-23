---
paths:
  - "R/scm.R"
  - "R/scm_fit.R"
  - "R/scm_rank.R"
  - "R/scm_frontier.R"
  - "src/scm.cpp"
  - "src/optim.cpp"
  - "src/inference.cpp"
---

# SCM（Abadie, Diamond & Hainmueller 2010; Ben-Michael, Feller & Rothstein 2022）

$$Y_{it}^{(0)} = \delta_t + \theta_t Z_i + \lambda_t \mu_i + \varepsilon_{it}$$

## 採用パターン・重み・共変量の設計不変条件

| コホート数 | コホート内ユニット数 | 経路 | 既定 |
|---|---|---|---|
| 1 | 1 | sharp | 通常の SCM |
| 1 | 複数 | sharp | **コホート内フルプーリング**（加重平均して SC 1 本） |
| 複数 | 任意 | staggered | `nu = "auto"` + `fixedeff = TRUE` |

- **同時採用の複数処置（BFR S.4.1）**: `Y_treat` は平均系列、生系列は `Y_treat_mat`（`T × n_treated`）。`n_treated > 1` の sharp fit では平均系列の低ノイズゆえ比率が過大になるため `mspe_ratio_pval()` を明示エラーにし、`scm_rank_pval()` へ誘導する
- **`scm_rank_pval()`（Abadie & L'Hour 2021 S.3.2）**: 置換時も採用時点をスロットに固定してイベントタイム構造を保つ。staggered は `pool_cohorts = FALSE` 必須
- **`augment_scm()`（BFR 2021, JASA）**: sharp SCM 限定（staggered は明示エラー）。リッジアウトカム回帰は無罰則切片（コントロール平均中心化）を持ち、符号契約は `tau_aug = tau_SCM - (m_tr_post - sum_j W_j * m_j_post)`（引き算、足し算禁止）。BFR 2022（staggered）とは別論文（`papers/Ben-Michael_et_al._2021` vs `_2022`）
- **ブロック粒度と質量 `m_g`（BFR A.2）**: `pool_cohorts` 既定 `TRUE`（参照実装 `augsynth` の `time_cohort = FALSE` と逆）で、pooled / `nu = "separate"` 両経路に効く。A.2 は和パラメータ化のため、`q_pool` に `m_g`、`q_sep` に `m_g^2` が掛かる（`.scm_block_target()`）。実装は平均スケールの系列を持ち、質量を係数側に出す
- **集計重み `weights`（Abadie 2021 S.8）**: `.scm_treated_weights()` で平均 1 に正規化し、`n_g` の代わりに `m_g` を使う。staggered の ①ブロックターゲット ②pooled 目的関数 ③`m_g * T_post` 集計 ④`scm_inference()` の per-unit 寄与 ⑤`scm_balance_frontier()` refit、および sharp のターゲット・`build_predictor_matrices(tr_weights=)` の全箇所で同じ重みを一貫して使う
- **`covariates` と `cov_method`**:
  - sharp SCM は `covariates` を明示エラーにする（ADH S.2.3 に従い `predictors = list(pred(...))` で balance する）
  - staggered `partial_out`（既定）は Clarke et al. (2024) S.2.2 の `.twfe_partial_out()`（`sdid.md` 参照）を借用
  - staggered `cov_method = "balance"`（BFR S.5.2 Eq. 11）: pooled ソルバ限定（`nu = "separate"` / `lambda_pen` / `donor_mspe_threshold` / `v_selection = "oos"` はエラー）。never-treated の平均・SD で標準化し `xi_sd` を掛けて `Xz'Xz` 項を足す（`xi` 既定は最終採用時点より前の never-treated 分散）。`xi` は balance 経路以外（sharp・`partial_out`）でエラーにする。**`augsynth` との 2 つの意図的差異**（論文準拠を採用）: ① `xi` 参照群が never-treated（augsynth は treated のラグアウトカム SD）② `q_*_rel` 正規化が脚注 7 の combined（augsynth はアウトカムのみ、そのため `q_sep_rel` が `nu = 0` で 1 にならず、アウトカム側の `q_sep`・`q_pool` が同時に下がることがある）
- **`nu = "auto"` の挙動**: 選択のない DGP ではノイズ平均化で `nu_hat ≈ 1/sqrt(J)` に寄る性質があるが、独自ヒューリスティックは入れず論文・`augsynth` 準拠を維持し、`scm_balance_frontier()` での確認に委ねる。`donor_mspe_threshold` / `lambda_pen` / `v_selection = "oos"` 併用時は既定が `"separate"` にフォールバックし原因引数を message で明示する（明示数値 `nu` はエラー）。staggered fit は解決後の `nu`・`qp_solver`・`fixedeff`・`pool_cohorts` を記録する
- **placebo 剪定と `mspe_threshold`**: `mspe_ratio_pval(mspe_prune =)` は事前 MSPE 比で placebo を落とす（`statistic = "post_mspe"` と組で使う）。`mspe_threshold`（比率分母の下限）とは別物。`plot.scm_placebo()` / `plot_data()` は既定 `mspe_prune = NULL` で検定のカットオフを引き継ぐ
- **sharp SCM の no-op 引数**: `control_group` と（`predictors = NULL` 時の）`scale_predictors` は既存共通呼び出し互換のため当面サイレント受理されている（未解決課題、`backlog.md`）

## 外側 V 最適化・内側 QP の数値契約

- **outcomes-only 自動ルーティング**: `predictors` が outcome のみの単年窓で `data` の全前処置期間をちょうど 1 回ずつ被覆する場合（`.is_outcome_only_spec()`）、`predictors = NULL` と同一経路に落とす（message あり、`v_optim = "multistart"` / `"bfgs"` 明示時は尊重）。被覆判定は `data` 内の前処置期間基準であり、単年ラグが不完全被覆のときは無言切替を防ぐため message を出す
- **multistart と `(V, W, loss)` の三つ組管理**: 非凸な外側最適化のため `v_optim = "auto"` は predictors 指定時 `"multistart"`、outcomes-only 時 `"coord_descent"`。評価子は KKT 検証済み解のみ許可し、座標降下が返す `(V, W)` の非整合を避けるため `(V, W, loss)` を必ず三つ組で追跡しベスト V から W を再導出しない。NNLS レスキュー（`use_nnls`）は multistart 内部限定で KKT ゲートは不変
- **k 次元因子形式 `(B, b)`**: 内側 QP は `B = diag(sqrt(V)) X0`（`k × N`）と `b = sqrt(V) * X1` で持ち、N×N Gram `Q = B'B` は作らない。`face_solve` は 3 経路（① `cheap_face` 用 bordered-KKT、特異時は諦めず厳密経路へ落とす ② `m > k` 用因子空間ヌル空間 solve ③ `m <= k` 用 Helmert）を持ち、**`wolfe_simplex()` から呼ぶときのみ `tight = true`（固有値 `1e-14`・Cholesky `1e-7`・特異値 `1e-11`）を使う**（`1e-9` 打ち切りだと `wolfe_simplex` の `1e-10` KKT 判定と衝突し自己再選択サイクルに入る。一方 `active_set_simplex` まで緩めると密な縮退面でピボット増になるため `wolfe_simplex` 限定）
- **outcomes-only 座標降下の Gram 形式**: `cheap_face && 2k >= N` では `active_set_simplex_gram()` を使い、Q のランク 1 更新＋暗黙表現でメモリ帯域詰まりを避ける
- **縮退 QP の非一意性**: `k < |support|` のとき最適解は面になり、`active_set` のタイブレークは浮動小数点丸めに依存するため **bit 再現を期待しない**。回帰検証は重み一致ではなく ① 自身の V での KKT 残差 ② 内側 QP 目的値 ③ 外側目的値で行う
- **`qp_solver = "auto"` とフォールバック順**: 予測子あり sharp と `2 * T_pre < N_co` の outcomes-only では `"wolfe"`、それ以外および R 側独自 QP 経路（`oos` / `lambda_pen` / `bfgs` / staggered、明示 `"wolfe"` はエラー）では `"active_set"`。`fit$qp_solver` に解決後の値を残し `mspe_ratio_pval()` が引き継ぐ。cold start は FISTA（support 絞り込み）→ active-set 仕上げの順（一様初期値からの active-set は O(N^3) 退化する）。`scm_inner_weights_cpp` も KKT 検証済み厳密解を返す。OSQP は外側座標降下の数万回呼出しに対し初期化コストが重いため不採用
- **リーク防止（`v_selection = "oos"` / `lambda_pen = "auto"`）**: outcomes-only では訓練期間（前半）のみで W を当て、後半検証期間の MSPE で V や `lambda_pen` を選び、最後に全期間/後半へ再フィットする（`v_window` は外側 V の評価行のみ制限し、`oos`・staggered と併用不可）

## staggered 推論・並列化

- **jackknife（`scm_inference()` 既定、BFR A.1）**: ν をフィット値に固定し全 N ユニットを LOO。高速化と整合のため ① ν=0 解（`pooling$separate_weights`）と最終解の両方で重み 0 のドナーは再推定を省略 ② 元の解から warm start ③ `solve_simplex_qp(lipschitz=)` に `lambda_max(FF')` を直接渡して毎回の `eig_sym` を省く。partial-out 済みアウトカム・balance の `Z` は固定
- **wild bootstrap（BFR S.5.3 / Otsu & Rai 2017 Eq. 8–9）**: 全ユニット（マッチしたドナーの残差も含む）に独立な黄金比 2 点乗数を掛ける。二方向（ユニット＋期間）残差から期間効果を除いて共通ショックの混入を防ぎ、**ATT で中心化するのは処置側の異質性項のみ**（全ユニット中心化すると SE が ATT に比例して膨らむ）。処置 3 未満で警告
- **OpenMP 並列（`src/inference.cpp`）**: placebo は（ドナー × multistart 候補）単位で並列化する。multistart 内部並列との二重並列は `omp_in_parallel()` でガードする
