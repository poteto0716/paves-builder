# 熱物性・引張弾性ワークフロー（PCFF / OpenMM）

繰り返し単位 SMILES と、1鎖および系全体の最小原子数を入力すると、条件を下回らない
最小の重合度（DP）と鎖数を決め、次を4つの独立seed（既定 `n=4`）で実行します。

1. PCFFで密度0.5 g/cm³の初期構造を作成
2. 800 KでNPT平衡化後、800→200 Kを40 K/ns（15 ns）でNPT冷却
3. **密度–温度曲線**の双曲線フィットからTgを決定
4. 密度フィットの低温・高温漸近勾配から線膨張係数を決定
5. 冷却中に保存した300 K構造を300 Kで再平衡化し、20 nsで2.0%まで一軸引張
6. 応力–ひずみの線形フィットから引張弾性率を決定

各物性はreplicaごとの値に加え、平均・標本標準偏差・標準誤差を出力します。

## 準備

```bash
conda env create -f workflows/thermomechanical/environment.yml
conda activate polypaves-thermomechanical
pip install dist/polypaves-0.3.0-cp313-cp313-linux_x86_64.whl
```

## 実行

```bash
cd workflows/thermomechanical

python thermomechanical.py new pmma \
  --monomer '*CC(C)(C(=O)OC)*' \
  --atoms-per-chain 1000 \
  --total-atoms 20000

python thermomechanical.py run projects/pmma
python thermomechanical.py status projects/pmma
```

### 共重合体をCLIから作る

主モノマーは `A`、追加モノマーは `--comonomer NAME=SMILES` で指定します。
`random` は組成を各鎖で厳密に合わせる `random-exact` の短縮名です。

```bash
# PMMA 70% / PS 30% ランダム共重合体
python thermomechanical.py new pmma_r_ps \
  --monomer '*CC(C)(C(=O)OC)*' \
  --comonomer 'B=*CC(c1ccccc1)*' \
  --architecture random \
  --fraction A=0.7 --fraction B=0.3 \
  --atoms-per-chain 1000 --total-atoms 20000

# PMMA 50% / PS 50% ジブロック共重合体
python thermomechanical.py new pmma_b_ps \
  --monomer '*CC(C)(C(=O)OC)*' \
  --comonomer 'B=*CC(c1ccccc1)*' \
  --architecture block \
  --fraction A=0.5 --fraction B=0.5 \
  --atoms-per-chain 1000 --total-atoms 20000

# 交互共重合体（fraction指定は不要）
python thermomechanical.py new pmma_alt_ps \
  --monomer '*CC(C)(C(=O)OC)*' \
  --comonomer 'B=*CC(c1ccccc1)*' \
  --architecture alternating \
  --atoms-per-chain 1000 --total-atoms 20000
```

確率的ランダム配列にする場合は `--architecture random-probabilistic` を使います。
各replicaでは `sequence_seed + replica index` を使うため、packingだけでなくランダム配列も
独立になります。確率配列でも全鎖が1鎖原子数の下限を満たせるよう、最小の繰り返し単位だけで
できた鎖を保守的な下限としてDPを決めます。

### PAVESファイルから作る

```bash
python thermomechanical.py new pmma_r_ps \
  --polypaves examples/pmma_random_ps.paves \
  --atoms-per-chain 1000 --total-atoms 20000

python thermomechanical.py new pmma_b_ps \
  --polypaves examples/pmma_block_ps.paves \
  --atoms-per-chain 1000 --total-atoms 20000
```

ファイル中の `monomer`、`terminator`、`architecture`、`block` などの化学指定を保持し、
`degree`、`chains`、密度、力場、温度、出力などはワークフロー設定で置き換えます。
`block A 70` / `block B 30` のような固定長ブロックは0.7 / 0.3の比率へ正規化し、
目標原子数を満たすDPへ自動調整します。単一の高分子・共重合体が対象で、混合物、溶媒、
外部構造、固定された明示配列はこの自動サイズワークフローでは受け付けません。

`--atoms-per-chain` と `--total-atoms` は水素を含む全原子数です。まずDP 1と2の
一鎖系をPAVESで実際に構築して、繰り返し単位を1個増やしたときの原子数を求めます。
その後、指定値以上となる最小DPと最小鎖数へ切り上げます。さらに、初期セル辺がPCFFの
非結合カットオフの2倍を超えない小系では、密度が1.5 g/cm³になってもOpenMMの
minimum-image条件を満たすまで鎖数を追加します（この保護密度は `project.json` の
`system.minimum_image_density_g_cm3` で変更可能）。実現した値は
`runs/01_build/sizing.json` に記録されます。

主なオプション:

| オプション | 意味 |
|---|---|
| `--terminator` | 末端基。既定 `*C`（メチル） |
| `--comonomer NAME=SMILES` | B以降のモノマー。複数回指定可能 |
| `--architecture` | `alternating`, `random`, `random-probabilistic`, `block` |
| `--fraction NAME=VALUE` | ランダム／ブロックの組成。合計1 |
| `--polypaves FILE` | PAVES化学指定ファイルから作成 |
| `--density` | 初期密度。既定0.5 g/cm³ |
| `--replicates` | 独立seed数。既定4 |
| `--forcefield` | PAVES `.ff`。既定 `examples/forcefields/pcff.ff` |
| `--platform` | `auto`, `CPU`, `CUDA`, `OpenCL`, `Reference` |
| `--test` | 全MD段階を数step～120 stepに短縮。接続確認専用 |

中断時は同じ `run` コマンドで未完了段階から再開できます。段階を明示してやり直すには
`--from cool`、単独実行には `--only analyze` を使います。

GPUを使用しない接続確認例:

```bash
python thermomechanical.py new smoke \
  --monomer '*CC*' --atoms-per-chain 100 --total-atoms 800 \
  --platform CPU --test
python thermomechanical.py run projects/smoke
```

`--test` は温度変化とひずみを極端に少ないstepへ押し込むため、得られるTg、膨張係数、
弾性率には物理的意味がありません。CLI、構築、OpenMM、状態分岐、フィット、集計が通るかだけを
確認するモードです。

## 計算条件

設定は作成した `project.json` にすべて保存されます。既定の本計算条件は次の通りです。

| 項目 | 既定値 |
|---|---:|
| replica | 4 |
| 力場 / MD | PCFF / OpenMM |
| 9-6分散の長距離補正 | 解析的tail補正あり |
| 初期密度 | 0.5 g/cm³ |
| タイムステップ | 1 fs（X–H SHAKE） |
| 圧力 | 1 atm |
| 800 K NPT事前平衡化 | 1 ns |
| 冷却 | 800→200 K、40 K/ns、15 ns |
| 冷却サンプリング | 10 psごと |
| Tgフィット用bin | 10 K |
| 弾性前300 K再平衡化 | **0.5 ns（50万step）** |
| 引張 | 工学ひずみ0→0.02、20 ns（ひずみ速度 1×10⁶ s⁻¹） |
| 弾性フィット範囲 | 0–2% |

引張軸は既定でxです。アモルファス系の方向平均も取りたい場合は `project.json` の
`elastic.axes` を `["x", "y", "z"]` に変更できます。この場合、各replicaの同じ300 K状態から
3方向へ独立に分岐し、方向平均をreplica値としてからreplica間平均を取ります。

## Tgと線膨張係数

10 Kごとに平均した密度 \(\rho(T)\) を次の滑らかな双曲線へフィットします。

\[
\rho(T)=\rho_0+m_g x+\frac{\Delta m}{2}
\left(x+\sqrt{x^2+w^2}\right),\quad x=T-T_g
\]

低温側の漸近勾配は \(m_g\)、高温側は \(m_g+\Delta m\) です。\(T_g\) は丸めた折れ曲がりの
中心、`transition_width_K` は丸みの幅です。等方的なアモルファス材料を仮定し、各領域の
線膨張係数を

\[
\alpha_L=-\frac{1}{3\rho(T_g)}\frac{d\rho}{dT}
\]

で求めます。低温側を `alpha_linear_glass_1_K`、高温側を
`alpha_linear_rubber_1_K` として出力します。

## 引張応力と弾性率

冷却ランプが300 Kへ到達した瞬間の `state_300K.xml` を全3方向1 atmの等方NPTで
**0.5 ns（50万step）再平衡化**し、その共通状態から各引張方向へ分岐します。引張軸のセルと座標を
連続的に伸ばし、横2方向はOpenMMのanisotropic barostatで1 atmに緩和します。軸応力は各
スナップショットでセルとfractional座標を微小に±変形し、PCFF全ポテンシャルエネルギーの
中心差分（PME項を含む）と理想運動項から評価します。0–2%の応力–工学ひずみを切片込みの
最小二乗直線で当てはめ、その傾きを引張弾性率とします。

## 出力

```text
runs/
├── 01_build/
│   ├── sizing.json
│   └── replica_01/ ... replica_04/   # PAVES/OpenMM初期系
├── 02_cool/
│   └── replica_XX/
│       ├── cooling.csv
│       ├── state_300K.xml
│       └── state_200K.xml
├── 03_elastic/
│   └── replica_XX/
│       ├── state_300K_equilibrated.xml
│       ├── stress_x.csv
│       └── state_x_2pct.xml
└── 04_analyze/
    ├── results.json
    ├── thermal_replicates.csv
    └── elastic_fits.csv
```

最終値は `runs/04_analyze/results.json` の `aggregate` にあります。`mean` が4 replica平均、
`std` が標本標準偏差、`sem` が標準誤差です。各フィットの `r_squared` と `rmse` も必ず確認し、
冷却速度・系サイズ・平衡化時間への収束性は別途評価してください。
