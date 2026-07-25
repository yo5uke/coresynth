# パフォーマンス特性

実測環境: Windows 11 x86-64 / R 4.6.1 / GCC 14.2.0 / Synth 1.1.10。両者とも outcomes-only 仕様（Synth は `special.predictors` に前処置全期間を指定）。再現スクリプト: `sandbox/24_bench_vs_synth_large.R`。

| シナリオ | 設定 | 単一フィット（中央値） | in-space プラセボ（全ドナー） |
|---|---|---|---|
| `large` | N_co=100・月次10年（T=120, T_pre=96, T_post=24） | Synth 81.4s / coresynth 0.03s → **~2,700x** | Synth 6,874s≈115分（5 ユニット実測→100 ユニットへ外挿）/ coresynth 0.28s → **~25,000x** |
| `small` | N_co=40・前処置2年+処置後5か月（T=29, T_pre=24, T_post=5） | Synth 6.76s / coresynth 0.0032s → **~2,100x** | Synth 250.4s≈4分（40 ユニット全数実測）/ coresynth 0.0095s → **~26,000x** |

両シナリオとも解は同一（ATT・重み cor=1.0000・pre-RMSPE が Synth と完全一致）— 同じ解に到達した上での速度差。要因は warm-started active-set QP・rank-1 Gram 更新 ＋ プラセボの OpenMP 並列（`scm.md` 参照）。倍率は N_co=40〜100 の範囲でほぼ一定（N_co が大きいほど優位性が縮むという傾向はない）。
