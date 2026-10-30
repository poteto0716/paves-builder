# 架橋エポキシネットワークの構築

溶媒を加える前の架橋モデルを、公開wheelのPAVES APIで生成します。
モノマーのSMILESからpackし、開始点を作成、加圧NPTと初期NPT緩和の後に
反応→最小化→NVT→NPTを繰り返します。反応後は局所Sageテンプレートで割り当てを更新します。
本体の実装はwheel内のコンパイル済みモジュールにあり、このスクリプトはAPIの実行例です。

リポジトリのルートから環境を準備します。

```bash
conda env create -f workflows/crosslinking/environment.yml
conda activate paves-reactive
pip install dist/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl
cd workflows/crosslinking
python crosslinking.py --project epoxy \
  --activation-fraction 0.1 \
  --max-reactions-per-cycle 50 --platform CUDA
```

既定モノマーは三官能TMP triglycidyl etherです。`--smiles`か`--smiles-file`で変更できます。
`--total-atoms`の既定は20000で、開始点作成時の追加Hも含めた目標です。
既定の`--atom-count-mode above`では、**20000原子を厳密に超える最小の整数分子数**へ
切り上げます。既定モノマー・開始点割合10%では421分子、初期化後20039原子です。
`nearest`（最も近い原子数）、`floor`（目標以下）、`exact`（厳密一致）も明示指定できます。
要求値と実現した原子数は`settings.json`と`metrics.json`に記録します。
packingではPAVESのnative局所緩和（`--pack-relax 50`、`--pack-postrelax 400`、PAVES既定と同じ）を行います。
以前は速度のため0にしていましたが、生の成長座標では嵩高い芳香族エポキシの結合が他分子の環を貫通し
（20000原子で約700本）、その後の最小化・MDでもほどけず、伸びた結合として最終構造まで残っていました。
緩和ありでは貫通・伸びた結合とも0本です（packingは数秒から約6分になります）。
packing後に環（3〜8員環）を貫通する結合が1本でもあれば停止します。
反応でも、新しい結合の線分が環（3〜8員環）を貫く組は反応させません。
反応後の最小化で断片が引きずられて環を通り抜けた場合（開環したエポキシに残るC-O結合など）は、
その反応をcycleから外してやり直します。外した数は`cycles.json`/`metrics.json`の
`ring_threading_withdrawn_reactions`、線分で除外した候補数は`ring_threading_rejected_candidates`です。
最終構造に環を貫通した結合が残っていれば、`final.*`を保存した後にエラーで停止します
（`metrics.json`の`ring_piercing_bonds_final`）。
以前は20000原子・開始点20%で1本残り、network_propertiesの最初のMDで結合長の検査に失敗していました。
packingの密度は既定0.5 g/cm³（`--density`）です。0.8で詰めると、嵩高い芳香族エポキシでは原子の重なりが最小化で取り切れず、
次の加圧NPTが4回中2回NaNで停止しました（0.5では0/4。加圧後は約1.1 g/cm³）。
packing直後のvoidを潰すため、初期NPTの前に300 K・200 MPa・50 ps（1 fs×50000 step）の
加圧NPTを行います。`--compression-pressure-mpa`と`--compression-steps`で変更でき、
`--compression-steps 0`で無効です。加圧後は座標・速度・boxをそのまま1 atmの初期NPTへ渡します。
初期NPTは300 K・1 atm・1 ns、各cycleは最小化、0.1 fs×10000 stepのNVT、
1 fs×1000 stepのNPTです。時間・温度・圧力・反応上限・距離はCLIで変更できます。
目標反応率は既定で99.5%（`--target-conversion 0.995`）、最大サイクル数は5000です。
目標に到達すれば早期終了し、未到達でも最大5000サイクルで最終構造を保存して終了します。
変更する場合は`--target-conversion`と`--cycles`を明示してください。

### 反応距離の上限と反応後の緩和

反応探索距離は5 Åから始まり、反応0件が2 cycle続くと1.2倍に広がり、最大7 Åまで伸びます（以前の既定は10 Å。7 Åを超えて作った結合が引き伸ばされたまま残ったため下げました）。
遠い対を無理に結合させると、伸びた結合や環の貫通などのおかしな構造ができやすくなります。
既定の上限は7 Åです。さらに抑える場合は`--max-reaction-cutoff 5`のように指定してください。
開始距離（`--reaction-cutoff`）を省略すると、開始距離は5 Åと上限の小さい方になります。
`--fixed-reaction-cutoff`を指定すると、距離は拡大しません。

各cycleの反応後の緩和は次の順で行います。

1. 最小化（`--minimization-iterations`、既定2000）。反応原子から`--reaction-minimization-radius`（既定12 Å）以内と
   更新した結合近傍の原子だけを動かし、それ以外は固定します。`0`で系全体を最小化します（以前の動作。系全体の熱的な構造も消えます）
2. 任意で`--post-reaction-md-steps`のNVT（0.1 fs刻み、反応のcommit前）
3. `--relaxation-stages`の各段

`--relaxation-stages`はJSON文字列、またはJSON/YAMLファイルで指定します。
指定すると`--nvt-*`／`--npt-*`の既定2段は使いません。
各段の書式は`ensemble`（NVT/NPT）、`steps`、`timestep`（fs）で、`temperature`（K）と
`pressure_atm`（NPTのみ）は省略できます。同じ内容の例を[relaxation.example.yaml](relaxation.example.yaml)に置いています:

```yaml
stages:
  - {ensemble: NVT, steps: 5000, timestep: 0.1}    # 0.5 ps: 新結合を小さい刻みでならす
  - {ensemble: NVT, steps: 10000, timestep: 0.5}   # 5 ps
  - {ensemble: NPT, steps: 20000, timestep: 1.0, pressure_atm: 1}  # 20 ps: 密度を戻す
```

```bash
python crosslinking.py --project careful --max-reaction-cutoff 5 \
  --post-reaction-md-steps 1000 --minimization-iterations 5000 --relaxation-stages relaxation.example.yaml
```

実際に使った段は`settings.json`の`cycle_md_stages`に、各cycleの結果は`cycles.json`の`md_stages`に残ります。

### 反応後の最終緩和

反応が終わったら（目標反応率、打ち切り、最大cycleのいずれでも）、最終構造を保存する前に
NVT 0.1 fs×1000 step（`--final-relaxation-steps`）→ NPT 1 atm・1 fs×1,000,000 step = 1 ns（`--final-md-steps`）を行います。
反応直後の構造は`reacted.structure.json`、緩和後を`final.*`として保存します。0で各段を省略できます。

### 反応が進まなくなったときの打ち切り

このモデルでは、反応できるのは活動状態をもつOH（開始点の数だけ）に限られます。
架橋が進んで活性点が網目に埋もれると、近くに未反応のエポキシがなくなり、反応0件のcycleが続きます。
探索距離が上限まで広がっても（固定距離ならその距離のまま）反応0件のcycleが
`--stop-after-idle-cycles`回（既定200）続いたら、その時点の構造を保存して終了します。
`0`を指定すると打ち切りません。終了理由は`metrics.json`の`stop_reason`に記録します
（`target_conversion`、`idle_at_maximum_cutoff`、`max_cycles`のいずれか）。

### 計算速度

高反応率で1 cycleが遅くなる主な原因と、その対策は次のとおりです。

- 新しい反応環境のテンプレート（局所フラグメントの電荷・パラメータ）計算：
  1 cycle内の未知テンプレートを`--template-workers`個のプロセスで並列計算します（既定4）。
  結果は逐次計算と同一です。1プロセスあたり約1 GBのメモリを使うので、メモリが少なければ減らしてください。
  計算済みテンプレートは`reaction_templates/`に保存して再利用します。
- 反応後の力場更新：反応点周辺の項だけを既存のOpenMM Systemの表から差し替えます。
  網目全体のOpenFF Topology/Interchangeは作り直しません。全体を作り直す場合と同じパラメータになることをテストで確認しています。
- 巨大分子のRDKit処理：検証と反応点探索で同じ変換結果を共有し、変化のない分子は再検査しません。
- 反応後の最小化：既定では反応点の周り（12 Å）だけを最小化します。系全体の最小化は熱揺らぎのため収束せず、
  毎cycle上限の反復数まで回っていました。

反応規則は同じディレクトリの`activation.yaml`と`reactions.yaml`です。
`--activation-fraction`／`--activation-count`を指定すると、`results/<project>/activation.yaml`の
`selection`と`fraction`／`count`の行もその値に書き換えて保存します（元のファイルは変更しません）。
これは中性OHへ活動状態を付けて伝播させる有効モデルで、帯電したカチオン機構そのものを
再現するモデルではありません。反応探索距離から実際の反応速度や障壁は求まりません。
ルールの定義は[YAML_GUIDE.md](YAML_GUIDE.md)を参照してください。

短いCPU動作確認:

```bash
OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python crosslinking.py --project smoke --total-atoms 500 --activation-count 1 \
  --platform CPU --compression-steps 2 --equilibration-steps 2 --cycles 1 \
  --nvt-steps 2 --npt-steps 2 --final-relaxation-steps 2 --final-md-steps 2 --max-reactions-per-cycle 1 \
  --reaction-cutoff 10 --max-reaction-cutoff 10 --minimization-iterations 100
```

`results/<project>/`に設定、反応履歴、cycleごとの集計、初期・最終構造を保存します。
容量を抑えるため、**OVITOトラジェクトリとcycleごとの構造ファイルは既定でOFF**です。
必要な場合だけ`--save-trajectory`を追加すると、初期・緩和後・各cycle終了時の
OVITO軌跡と、各cycleの構造を保存します。
`final.mol2`には部分電荷と結合・結合次数、`final.structure.json`にはさらに
安定原子ID、活動状態、周期セルを保存します。吸水計算にはJSONを渡します。
周期境界をまたぐ架橋を通常の有限分子用MOL2 importで読み直さないでください。

同名projectがある場合は停止します。明示的な`--force`はそのprojectの結果を削除します。
架橋構築の途中再開には対応していません。反応履歴と集計は軌跡OFFでも更新します。
CPUの短い実行はパイプラインの確認で、目標架橋率や樹脂の平衡を保証しません。

作成した構造は[架橋ネットワークの物性評価](../network_properties/README.md)へ渡せます。
800→100 Kの冷却とTg/CTE、300 Kでの熱・弾性・熱輸送物性、任意の中性溶媒の
挿入自由エネルギーと拡散係数を評価します。
溶媒をpackしてMDだけ行う場合は[溶媒packと吸水MD](../sorption/README.md)を参照してください。
