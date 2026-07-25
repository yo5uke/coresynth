---
paths:
  - "R/scm.R"
  - "R/scm_fit.R"
  - "src/scm.cpp"
  - "src/optim.cpp"
---

# SCM（Abadie, Diamond & Hainmueller 2010）

$$Y_{it}^{(0)} = \delta_t + \theta_t Z_i + \lambda_t \mu_i + \varepsilon_{it}$$

重み推定は入れ子最適化: 外側 = V（予測子重み）の座標降下 / L-BFGS-B / multistart、内側 = simplex QP（FISTA + active-set）。推論は `mspe_ratio_pval()`（sharp、MSPE 比率順列検定）と `scm_inference()`（staggered 専用 wild bootstrap）に分離。検証ツール: `placebo_in_time()`（in-time placebo backdating）・`loo_donors()`（V 固定ドナー LOO）。

**Staggered**: `predictors = NULL` 限定（`covariates` の partial-out のみ対応）。理論根拠は Ben-Michael, Feller & Rothstein (2022, JRSS-B)。原著 SCM は単一処置ユニット前提のため、staggered 拡張はこの論文に依拠する — コホート内平均（完全プーリング）は同一採用時点内で DGP 異質性項が消えることから正当化される。コホート間は `nu = NULL`（既定、per-cohort 独立 V 最適化）または `nu ∈ [0,1]` / `"auto"` の **partially pooled SCM**（コホート横断の pooled 目的と各コホート separate 目的を ν で線形結合、全コホートの重みを同時最適化）。`fixedeff = TRUE` で intercept shift（各ユニットを自身のプリ平均でデミーンした weighted DiD 形式；報告される `Y_synth` はこの shift を戻した生スケール）。

## 内側 QP・外側最適化

- **内側 simplex QP**: FISTA に適応リスタート（O'Donoghue & Candès 2015）を導入している。事前アウトカムの共線性で Q 行列が悪条件になりやすく、リスタートなしだと反復数が条件数に比例して膨張するため。解は quadprog と高精度で一致する
- **外側 V 最適化は非凸**。単一スタート（座標降下・BFGS）は予測子つきの経路で悪い局所解に停留しうる（ADH 原論文の重みを再現できないケースを確認済み）。既定 `v_optim = "auto"` は predictors 指定時 `"multistart"`、outcomes-only 時 `"coord_descent"`（outcomes-only では局所解問題が実害として出ていないため軽い方を使う）
- **outcomes-only 標準形の自動ルーティング**: predictors が「outcome 変数のみ・各 `pred()` が単年・前処置全期間をちょうど 1 回ずつ被覆」（`.is_outcome_only_spec()`、単年窓では op 不問・順序不問）なら X0/X1 が前処置アウトカム行そのものになるため、`predictors = NULL` と同一の outcomes-only 経路に落とす（message あり）。Synth はこの per-year special.predictors 指定を outcomes-only の標準形とするため、同じモデルの 2 通りの書き方が同一の解を返すことを保証する措置。`v_optim = "multistart"` / `"bfgs"` の明示指定時はルーティングせず予測子経路を尊重する
- **multistart の設計**: 決定的な複数スタート（uniform / one-hot / 固定シード Dirichlet 乱数）を warm-start QP でスクリーニングし、上位候補のみ Nelder-Mead↔座標降下で磨く。**評価子は KKT 検証済みの解のみ許可**する — 未収束の内側解を通すと外側目的が不正確になり最適化が迷走するため。座標降下が返す (V, W) は暗黙に非整合ペア（W はスイープ中の中間解）なので、**(V, W, loss) は必ず三つ組で追跡**し、ベスト V から W を再導出しない。placebo は（ドナー×候補）をまとめて並列化する — ドナー単位で並列化すると最も遅いドナーが1コアを占有し他コアが遊ぶため。OpenMP のネスト（multistart 内部の並列 × placebo の並列）は `omp_in_parallel()` でガードし二重並列によるオーバーサブスクリプションを避ける
- **NNLS レスキュー**（`nnls_rescue` フラグ、既定 false・multistart 内部限定）: 縮退したドナー幾何（ほぼ重複する予測子行）では厳密な active-set 解法が毎回コールドスタート FISTA にフォールバックし遅い。Lawson-Hanson NNLS で非ゼロサポートだけを先に推定してから同じ厳密ソルバに渡すことで解消する。**KKT 検証ゲートは変わらない**ので厳密性は落ちない。plain `coord_descent`・oos・staggered・placebo の通常経路ではこのフラグは常に false で、経路自体は変えていない
- **面 QP ソルバは 2 経路ある**（`cheap_face` フラグで切替）: 予測子指定 / multistart の経路は**和制約のヌル空間（Helmert 基底）で解く**必要がある — multistart は one-hot V 由来のランク落ち G を頻繁に踏むため、bordered-KKT 直接 solve だと失敗しやすい。加えて分岐基準をデータスケールに固定しないと `scale_predictors` のスケール不変性が壊れる。一方 outcomes-only はこのリスクがないため、**より軽い bordered-KKT 直接 solve**（`cheap_face = TRUE`）を使い壁時間を節約している。**この使い分けを逆にしない** — 予測子経路に cheap_face を使うとスケール不変性が壊れ、outcomes-only にヌル空間ソルバを使うと不要に遅い
- **`v_window`**（sharp SCM 限定）: 外側 V フィットの評価行だけを制限する。X（予測子行列 or outcomes-only の全プリ行）と報告される `loss`（全プリ期間）は不変。`v_selection="oos"` と staggered とは併用不可
- **OOS V 選択**（`v_selection="oos"`, outcomes-only, Abadie 2021 §3.2 準拠）: プリ期間を訓練半分と検証半分に分割し、W(V) は**訓練半分のみ**でフィットする（検証行を QP に含めると V が検証行に質量を寄せて検証 MSPE を人為的にゼロにできるリークが起こる）。V は検証 MSPE 最小化で選び、最終 W は V* を用いて直近の半分に再フィットする
- **予測子 SD スケーリング**（`scale_predictors = TRUE` 既定）: 各予測子を全ユニット横断 SD で除算する。スケールの大きい予測子が V 加重損失を支配するのを防ぐため。W 自体はスケール不変

## staggered・推論

- **partially pooled の解法**: コホート単位の巡回ブロック座標降下（他コホート固定下では標準 simplex QP に帰着）。イベントタイム整列時の欠損ラグはゼロ詰め。`nu` は `donor_mspe_threshold`/`lambda_pen`/`v_selection="oos"` と併用不可（pooled 経路は outcomes-only・uniform lag weights を前提とするため）
- **`control_group="clean"`（既定）の下方バイアス**: future adopter をドナーに含めるため、後期コホートの処置後アウトカムが前期コホートの反実仮想を汚染し ATT が下方バイアスを持つ。`"never_treated"` にすると消えるが、ドナープールが縮む。既定は後方互換のため変更していない
- **wild bootstrap**（`scm_inference()`、BFR S.5.3 / Otsu & Rai 2017）: ATT を per-treated-unit 寄与の加重平均に書き直し、黄金比 2 点乗数で帰無分布を生成する。CI・p 値はこの経験分布から求める（正規近似ではない）。`fixedeff` fit では寄与を DiD 形式（プリ平均でデミーン）で計算する。処置ユニットが 5 未満だと警告する
- **OpenMP 並列化**（`src/inference.cpp` の `scm_placebo_cpp`/`scm_placebo_x_cpp`）: per-control-unit の placebo イテレーションは書き込み先が独立インデックスのため完全並列化できる。壁時間は最も遅い1ドナーで律速される（コア数を増やしても短縮しない）
- **OSQP は不採用**: 座標降下が simplex QP を1フィットあたり数万回呼ぶ構造上、OSQP の初期化・前処理コストが累積して不利。FISTA + simplex 射影の方が N_co ≤ 50 で優位
