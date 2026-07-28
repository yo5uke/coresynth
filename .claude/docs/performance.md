# パフォーマンス特性

実測環境: Windows 11 x86-64 / R 4.6.1 / GCC 14.2.0 / Synth 1.1.10。両者とも outcomes-only 仕様（Synth は `special.predictors` に前処置全期間を指定）。再現スクリプト: `sandbox/24_bench_vs_synth_large.R`。

| シナリオ | 設定 | 単一フィット（中央値） | in-space プラセボ（全ドナー） |
|---|---|---|---|
| `large` | N_co=100・月次10年（T=120, T_pre=96, T_post=24） | Synth 81.4s / coresynth 0.03s → **~2,700x** | Synth 6,874s≈115分（5 ユニット実測→100 ユニットへ外挿）/ coresynth 0.28s → **~25,000x** |
| `small` | N_co=40・前処置2年+処置後5か月（T=29, T_pre=24, T_post=5） | Synth 6.76s / coresynth 0.0032s → **~2,100x** | Synth 250.4s≈4分（40 ユニット全数実測）/ coresynth 0.0095s → **~26,000x** |

両シナリオとも解は同一（ATT・重み cor=1.0000・pre-RMSPE が Synth と完全一致）— 同じ解に到達した上での速度差。要因は warm-started active-set QP・k 次元因子形式の内側 QP ＋ プラセボの OpenMP 並列（`scm.md` 参照）。倍率は N_co=40〜100 の範囲でほぼ一定（N_co が大きいほど優位性が縮むという傾向はない）。

## 予測子経路（`predictors` 指定）

上表は outcomes-only。予測子指定の経路は `Q = X0'VX0` の rank が高々 k のため**予測子が少なくドナーが多いほど内側 QP が縮退**し、旧実装（N×N Gram + O(m³) face solve）ではここが支配的だった。因子形式への置き換え後の実測（N_co=240・T_pre=20・T_post=10・16 コア、再現: `sandbox/28_tier1_qp_refactor_check.R` と同じ DGP）:

| ケース | 単一フィット | in-space プラセボ（240 ドナー） |
|---|---|---|
| outcomes-only | 0.07s → 0.07s | 1.89s → **1.14s** |
| k=1（`pred(y, 1:20)`）・既定 multistart | 3.78s → **<0.01s** | 315s → **0.03s** |
| k=3・既定 multistart | 5.79s → **0.05s** | 463s → **2.11s** |
| k=3・`coord_descent` | 0.47s → **0.02s** | 75.3s → **0.37s** |

**k が小さいほど遅い**という反直感的な性質が旧実装の特徴だった（rank 落ちが深いほど active set が全ドナーから始まり O(m³) の face を踏み続けるため）。現在は k=1 が最速。
