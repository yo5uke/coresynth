---
paths:
  - "R/gsc.R"
  - "src/gsc.cpp"
---

# GSC（Xu 2017）

$$Y_{it} = \delta_{it} D_{it} + x_{it}'\beta + \lambda_i' f_t + \varepsilon_{it}$$

## 推定と加法的固定効果 `force`

- **集中化 `beta` 更新（Bai 2009）**: 共変量あり反復では素朴な ALS ではなく集中化形 `beta = (sum X' M_F X)^(-1) sum X' M_F Y` を使う（共変量と因子の相関下での収束確保、既定 `max_iter = 1000`, `tol = 1e-8`）。未収束時は `.gsc_warn_unconverged()` で警告
- **Step 1 / Step 2 分離**: `.gsc_fit_controls()` が対照 IFE（`gsc_ife_cpp` に 0 列処置渡し）、`.gsc_counterfactual()` が処置側（非 export `gsc_treated_cpp`）。staggered `"never_treated"` では対照が全コホート共通のため IFE を 1 回だけ当てて全コホートのローディングを算出する（`"clean"` はコホート別）。staggered fit は `control_group` を保持する
- **`force`（Xu 2017 Remark 3）**: `"none"`（既定）/ `"unit"` / `"time"` / `"two-way"`。対照を within 変換してから Bai 反復を回し、残差から `mu`・`alpha_co`・`xi` を復元し、処置側は事前期間のみで `Y - X*beta - mu - xi_t` を `[1, F_tilde]`（または `F_tilde`）に回帰して `(alpha_i, lambda_i)` を得る（処置を全期間でデミーンすると事後期間が混ざるため禁止）。`force != "none"` では `r = 0` も許可・CV 候補に含める。**`gsynth` との差**: `"none"` のみ `gsynth` が総平均 `mu` を推定する分だけ異なり、coresynth は Xu 式(1) どおり `mu` を持たない
- **因子の正規化と CV（Algorithm 1）**: `F = U_k * diag(D_k)` は `F'F = T*I_r` を満たさないが、射影ベースの反事実には影響しない。CV（`.gsc_cv_r()`）は MSPE(r) の分母を `T_pre` のみにし（`force = "none"` では候補 `r >= 1`）、staggered では全処置が未処置のグローバル事前窓で 1 回だけ選択する

## 推論契約

- **`gsc_boot()`（sharp 限定・Xu 2017 Algorithm 2）**: 対照 LOO 予測誤差ベクトルを処置側に、対照 in-sample 残差ベクトルを対照側に**時系列ごと**再抽出する（系列相関保持）。CI は basic 法 `ATT_obs - q(1 - alpha/2)` 〜 `ATT_obs - q(alpha/2)`（複製は帰無で生成される推定誤差分布のため percentile 字面は採らない）、p 値は `(1 + #{|ATT*| >= |ATT|}) / (B + 1)`
- **`gsc_inference()`（`.unit_inference()` 共通ドライバ）**: 処置も層別抽出し、staggered では処置をプールしてコホートを再構築、対照は 1 回だけ引いて全コホートで共有する。jackknife で処置を落とす回は `ctx$cache` の対照 IFE を再利用する。`N_tr = 1` はエラー、staggered jackknife は全コホート `N_tr,g >= 2` 必須
- **`conformal_inference()`**: 共変量つき fit では `beta` と処置平均共変量を含めて再推定する（ドロップ禁止）
