# パフォーマンス特性

実測環境: Windows 11 x86-64 / R 4.6.1 / GCC 14.2.0 / Synth 1.1.10。両者とも outcomes-only 仕様（Synth は `special.predictors` に前処置全期間を指定）。再現スクリプト: `sandbox/24_bench_vs_synth_large.R`。

| シナリオ | 設定 | 単一フィット（中央値） | in-space プラセボ（全ドナー） |
|---|---|---|---|
| `large` | N_co=100・月次10年（T=120, T_pre=96, T_post=24） | Synth 69.9s / coresynth 0.020s → **~3,500x** | Synth 5,817s≈97分（5 ユニット実測→100 ユニットへ外挿）/ coresynth 0.39s → **~15,000x** |
| `small` | N_co=40・前処置2年+処置後5か月（T=29, T_pre=24, T_post=5） | Synth 5.78s / coresynth 0.0028s → **~2,100x** | Synth 209.6s≈3.5分（40 ユニット全数実測）/ coresynth 0.012s → **~17,500x** |

**計測ルール**: 倍率はマシン状態で大きく揺れる（50% 程度の差が出うる）ため、版間の速度比較は**必ず同じ実行内で交互に計測し最小値で比較する**。性能に触る変更の後は `sandbox/24_bench_vs_synth_large.R` の coresynth 側だけでも回し、上表と桁が合うか確認する（テストは数値結果しか見ないため速度後退を検出できない）。

両シナリオとも解は Synth と一致（ATT・重み cor=1.0000・pre-RMSPE 完全一致）。要因は warm-started active-set QP・k 次元因子形式の内側 QP ＋ プラセボの OpenMP 並列（`scm` ルール参照）。

## 予測子経路（`predictors` 指定）

k 次元因子形式の実装による改善実績（N_co=240・T_pre=20・T_post=10・16 コア、再現: `sandbox/28_tier1_qp_refactor_check.R` と同じ DGP）:

| ケース | 単一フィット | in-space プラセボ（240 ドナー） |
|---|---|---|
| outcomes-only | 0.07s → 0.07s | 1.89s → **1.14s** |
| k=1（`pred(y, 1:20)`）・既定 multistart | 3.78s → **<0.01s** | 315s → **0.03s** |
| k=3・既定 multistart | 5.79s → **0.05s** | 463s → **2.11s** |
| k=3・`coord_descent` | 0.47s → **0.02s** | 75.3s → **0.37s** |

## その他手法の速度基準

- **SCM staggered jackknife**: N=200 で約 3.0s（zero-weight ドナー除外＋warm-start＋`lipschitz` 直接渡しによる）、N=400 で約 80s（1 回の pooled fit 0.38s、FISTA の O(N²) 勾配が律速）。
- **SDID（`sdid_inference(method = "placebo")` 対 synthdid 0.0.9）**: T=60・T_pre=45 のパネル（`sandbox/32_sdid_placebo_vs_synthdid.R`）で 1 placebo あたり coresynth 0.0010–0.0011s（synthdid の 250–930 倍高速）。N_co=300 でも FISTA の横長回避により 4–8 倍高速（N_co=299・N_tr=1 で 1.3s、N_tr=5・200 draws で 0.6s）。
- **TASC**: N=200・T=50 で約 0.2–0.34s（push-through 恒等式＋対角 R ベクトル持ち回り、中央値 350–400 EM 反復）。
