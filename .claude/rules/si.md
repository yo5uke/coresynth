---
paths:
  - "R/si.R"
  - "src/si.cpp"
---

# SI（Agarwal et al. 2025）

$$Y_{ti}^{(d)} = \sum_\ell u_{t\ell} v_{i\ell} \lambda_{d\ell} + \varepsilon_{ti}^{(d)}$$

SI-PCR: $Y_{pre,co}^\top$ の truncated SVD → 疑似逆行列でドナー重みを計算する。k の既定は Gavish–Donoho の最適ハード閾値（下記）。推論は `si_inference(method="bootstrap"/"jackknife"/"jackknife_global")` — sharp・staggered・multi-arm・staggered×multi-arm すべてに対応。

**既知の制限**: 共変量は未対応（sharp/staggered とも）。原著は本体では共変量なしだが、**Appendix E に時不変の共変量を事前アウトカムの下に積んで重みを学習する拡張がある**（`Z_i = [Y_pre,i; X_i]`、Theorem 2 は T0 を T0+K に置き換えて成立）。以前ここに「原著が covariate なしの定式化」と書いていたのは不正確だった。実装する場合はこの Appendix E 形式（ADH の予測子と同型）になる

**推論はユニット再抽出**（1.1.0 で変更、`gsc.md` の同名節と共通のドライバ `.unit_inference()`）: 旧実装は処置を固定して対照だけを引いており、層別にすると被覆が N_tr=5 で 0.43→0.88、20 で 0.73→0.93、staggered（3 コホート×5）で 0.47→0.84 に戻る。Agarwal et al. は bootstrap を提案しておらず（Appendix G の漸近正規性は donor 部分集合を選ぶ修正版推定量が前提で実用的でない）、ユニット再抽出は fect・Clarke Algorithm 2 と同じ規約を採った。**N_tr=1 は明示エラー**（sharp なら `conformal_inference()` へ誘導）、staggered の jackknife はどのコホートにも処置 2 以上が必要（Clarke Algorithm 3）。multi-arm は**アーム内で**処置を層別に引き（アーム間はプールしない）、どれかのアームの処置が 1 つならエラー。各再推定は fit が使った k（コホート別）をそのまま使う。`.si_unit_estimate()` はコホート内の全アームの処置列を 1 回の `si_pcr_cpp` で解く（W は処置列に線形なのでアーム別に解くのと同一）

**staggered×multi-arm の per-cohort jackknife の旧バグ**: コホート内アームの**加重和**（`Σ_d w_gd τ_gd`）の分散に `(w_g/W)²` を掛けており、SE が w_g（= n × T_post、数十）倍に膨らんでいた。共通ドライバではコホートの加重**平均**の分散を使う

**Multi-arm (K>1)**: ユニット因子 $v_{i\ell}$ が arm に依存しないという性質を使い、コントロール arm (d=0) の SVD 基底を全 K 処置 arm で共有する。`d` が多値（`max(d) > 1`）なら自動的に multi-arm 経路へ dispatch する。

**Staggered 設計**: `Y_pre_co_g`（T_pre_g × N_co_g）のみで SVD を取る — future-treated も T_pre_g 以前は未処置なので汚染がない。**Staggered × multi-arm** はコホート g × アーム d の (g,d) セルごとに SI-PCR を行い、コホート内で SVD 基底を共有する。集計式:

$$\hat\tau = \frac{\sum_g\sum_d N_{tr,g,d} T_{post,g} \hat\tau_{g,d}}{\sum_g\sum_d N_{tr,g,d} T_{post,g}}$$

**Multi-arm 推論の注意**: bootstrap はコントロール列を全 arm・全コホートで**共有サンプリング**する — SVD 基底が shared donor pool を前提にしているため、arm ごとに独立リサンプルすると基底の整合性が崩れる（以前は staggered でコホートごとに独立に引いていたが、never-treated 既定では全コホートが同じ対照を使うので 1 回に統一した）。sharp multi-arm には `jackknife_global` が無い（staggered 専用のため）。

**k の選択**: Agarwal et al. は CV・universal SVHT・elbow の 3 つをデータ駆動の方法として挙げる。採用したのは 2 番目 — `tau = omega(beta) * median(singular values)`、`beta = min(dim)/max(dim)`、`omega` は原論文 Eq.5 の閉形式近似（Gavish & Donoho 2014 は 0.001 < β ≤ 1 で真の係数との差 ≤ 0.02 と報告。2026-09-23 に原典で係数・手順とも照合済み。以前ここにあった「誤差 0.0001」は誤記だった）。**SVD 1 回で済む**のが決め手（CV だと候補ごとに PCR フィットが要る）。旧既定 `floor(sqrt(min(T_pre, N_co)))` は上記いずれでもない独自式だった。既知ランクの DGP で r=1〜4 を完全復元（旧式は常に 4 を返していた）。staggered/multi-arm はコホートの `Y_pre_co_g` に対して同じ閾値を当て、`min(T_pre_g, N_co_g) - 1` で上限を切る。**したがって staggered には単一の k が存在しない** — GSC の `r` と MC の `lambda` がパネル全体で 1 回の探索なのと違う。`fit$k` は採用時点を名前に持つコホート別の整数ベクトル（`.si_k_by_cohort()`、multi-arm はコホート内で SVD 基底を共有するので arm 間で同値）。sharp はスカラーのまま。v1.0.0 の途中までは**引数（既定 NULL）をそのまま格納**していて選んだ値を報告していなかった
