# 単分子から架橋膜までの物性評価

2つのコマンドに分けて計算します。

1. **単分子プロジェクト**（`monomer`）: モノマーをpackし、温度T（既定363.15 K = 90 ℃、`--temperature`）で緩和します。
   続けて熱履歴（800→100 K）でTg、300 Kの物性、最後にTの物性（環の数密度・配位数・MSDなどを含む）を求めます。
   Tの最後のNPT構造`replica_XX/<T>K/relaxed_<T>K.structure.json`が、架橋計算の出発構造です。機械特性は計算しません。
2. **架橋プロジェクト**（`crosslink`）: 単分子プロジェクトのディレクトリと開始点の密度を指定します。
   密度ごとに`<単分子>/crosslinks/<名前>/`を作り、架橋レプリカNは**単分子レプリカNのT緩和構造**に開始点を作って、
   同じ温度Tで架橋します。その後、熱履歴（Tg）と300 Kの物性を既定で求めます。
   機械特性（弾性）は時間がかかるため、膜ができた後に`mechanical`コマンドで任意に実行します。

保存済みの乾燥ネットワークを評価する場合は`network`コマンドを使います（Tg・300 K・溶媒）。
架橋後の結合、部分電荷、安定原子IDを保持し、その後は固定トポロジーのSage/OpenMMで計算します。
既存のPCFF鎖構築用`thermomechanical.py`とは別の入力経路です。

計算する物性の一覧と、それぞれの意味・求め方は[PROPERTIES.md](PROPERTIES.md)にまとめています。

## 実行

架橋構築と同じ[Python 3.11 / OpenFF環境](../crosslinking/environment.yml)と
公開wheel `dist/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl`を使用します。
リポジトリのルートで実行してください。

```bash
conda activate paves-reactive
NP=workflows/network_properties/network_properties.py

# 1. 単分子：4レプリカ（独立なpacking）、90 ℃
python $NP monomer workflows/network_properties/output/TMPTGE --replicates 4 --seed 31415 --platform CUDA
python $NP submit workflows/network_properties/output/TMPTGE --monomer-only   # 単分子だけ（または run）
python $NP status workflows/network_properties/output/TMPTGE                   # 全レプリカの完了を確認

# 2. 架橋密度を並べて指定すると、密度ごとに子プロジェクトを作る（単分子の計算前に作ってもよい）
python $NP crosslink workflows/network_properties/output/TMPTGE \
  --activation-fraction 0.02 0.05 0.1 --max-reactions-per-cycle 50 --solvent water=O
python $NP submit workflows/network_properties/output/TMPTGE          # 単分子と全架橋をまとめて投入

# 3. 任意：機械特性（架橋膜の300 K構造から）
python $NP mechanical workflows/network_properties/output/TMPTGE/crosslinks/af0.05
python $NP submit workflows/network_properties/output/TMPTGE --mechanical   # GPUへまとめて投入

python $NP status workflows/network_properties/output/TMPTGE   # 単分子と各架橋の進捗をツリー表示
python $NP report workflows/network_properties/output/TMPTGE   # レポートと架橋シリーズを書き直す（MDなし）
```

ディレクトリ構成:

```text
TMPTGE/                               単分子プロジェクト
  project.json  definitions/          開始点・反応ルール（activation.yaml / reactions.yaml）
  replica_01/
    build/                            pack → 200 MPa加圧 → T・1 atmで1 ns緩和（開始点なし、反応なし）
    00_minimize … 06_production_300K/ 熱履歴（800→100→300 K）と300 K NPT
    07_nvt_300K/
    363.15K/                          300 K → T → 1 ns緩和 → 2 ns NPT（物性・環の統計）→ 2 ns NVT
      relaxed_363.15K.structure.json  ★架橋の出発構造
      03_production/reactive_rings/   環の数密度・配位数・クラスター・MSD
  report.md  property_summary.csv
  crosslink_series.csv / .png         架橋密度ごとの比較（Tg、反応率、密度、残膜率、網目密度…）
  crosslinks/
    af0.05/                           架橋プロジェクト（開始点の割合0.05）
      project.json  definitions/
      replica_01/                     ← 単分子replica_01の★から開始
        crosslinking/results/network/ 開始点作成 → opt → NVT → NPT → 反応cycle → 最終緩和
        network_origin.json           出発構造（単分子レプリカ・ハッシュ）、反応率、開始点密度
        00_minimize … 07_nvt_300K/    熱履歴（Tg）と300 K物性
        mechanical/                   任意（mechanicalコマンド）
      report.md
    n8/                               --activation-count 8 の場合
```

### 単分子プロジェクト（`monomer`）

| 段階 | 既定 |
|---|---|
| 構築 | packing（0.5 g/cm³、`--packing-density`）→ T・200 MPa・50 psの加圧NPT（`--compression-pressure-mpa`、`--compression-steps`）→ T・1 atm・1 ns（`--equilibration-steps`） |
| 熱履歴 | 300→800 K加熱、800 K 1 ns、800→100 K冷却（40 K/ns）でTg、100→300 K |
| 300 K | NPT 10 ns（密度・Cp・CTE・圧縮率…、環の数密度）、NVT 10 ns（直接Cv・誘電） |
| T | 300 K NPTの続きからTへ20 K以下の刻みで昇温（各100 ps）→ 1 ns緩和 → NPT 2 ns → NVT 1 ns緩和 + 2 ns |

Tの物性は`replica_XX/<T>K/result.json`と`report.md`の「<T> K monomer liquid」に出ます。
環（エポキシ・オキセタンなど、開始点ルールの`target_smarts`に一致する未反応環）について、
300 KとTの**単位体積あたりの環の数**（nm⁻³とmmol/cm³。環の数 ÷ その温度のNPT平均体積）を
それぞれ「300 K thermodynamics」と「<T> K」の表に加えます。
TのNPT観測（1 psごと）からは、環の配位数・接触グラフの最大クラスター分率（4 Åと6 Å、`ring_contact_cutoffs_A`）、
環の重原子重心のMSD・拡散係数・MSD指数αも求めます（定義は[crosslinking/README.md](../crosslinking/README.md)）。
拡散係数は、MSD指数αが0.8〜1.2（拡散域）の場合だけチェック通過として集計します。
`project.json`の`t_equilibration_ps`、`t_production_ps`、`t_nvt_production_ps`で時間を変えられます。

既定のモノマーはTMP triglycidyl ether（`--smiles`で変更）です。`--chemistry oxetane`で
オキセタン用のYAML（[crosslinking/oxetane/](../crosslinking/oxetane/)）と既定モノマーOXT-221
（`CCC1(COCC2(CC)COC2)COC1`）を使います。ほかの環や反応は`--activation-yaml`／`--reactions-yaml`で渡します。
ルールは`definitions/`へコピーし、SHA-256で変更を検査します。表の見出しは環に合わせて
「Epoxide …」「Oxetane …」になります。

単分子の液体は800 K・1 atmで沸騰・空洞化しうるので、各段階の終了時密度を
`dry_result.json`の`thermal_history_density`に記録し、800 Kの密度が開始時の0.7倍を下回ると
Tg/CTEと300 K NPTのチェックを不合格にします（後述）。必要なら`--high-temperature`で下げられます。

### 架橋プロジェクト（`crosslink`）

```bash
python $NP crosslink <単分子dir> --activation-fraction 0.05 0.1     # 子プロジェクト af0.05, af0.1
python $NP crosslink <単分子dir> --activation-count 4 8 --name ...  # n4, n8（--nameは1密度のときだけ）
```

- 架橋レプリカNは単分子レプリカNの★構造から始めます（番号で固定）。レプリカ数の既定と上限は単分子のレプリカ数です。
  単分子レプリカNが未完了なら、`run`はエラーで止まり、`submit`は完了まで待ちます（単分子が失敗した番号は`blocked`）。
- 反応温度は単分子のTに固定です（変えたいときは別の単分子プロジェクトを作ります）。
- 開始点を作った後、opt → NVT 0.1 fs × 10000 step → T・1 atmのNPT 100 ps（`--equilibration-steps`）で緩和してから
  反応cycleへ進みます。開始点で足すHの分だけ原子数は増えます。
- seedの既定は単分子と同じなので、架橋レプリカNのseedは単分子レプリカNと同じです。開始点の選択順が密度間で共通になり、
  割合を増やすと開始点の集合が入れ子になります（密度間の比較の雑音が減ります）。`--seed`で変えられます。
- 単分子の`project.json`を作成後に変えると、架橋側は実行を拒否します。単分子レプリカの★構造が架橋後に変わった場合も
  （`reseed`など）、`network_origin.json`のハッシュと一致しないので拒否します。
- 残膜率`structure.film_retention_volume_ratio`は、架橋膜の熱履歴後300 K NPT平均体積 ÷ **同じ番号の単分子レプリカの
  300 K NPT平均体積**です。両方とも300 Kなので熱膨張は入りません（`--remove-sol-fraction`なら分子はゾル削除後のゲル）。
- 開始点数・割合（実際の値）、架橋開始直前の体積あたりの開始点数、反応率、活性末端1個あたりの連鎖長などを
  `network_origin.json`と「Network construction」に出します。

反応の設定は従来どおり指定できます: `--cycles`（既定5000）、`--target-conversion`（既定0.995）、
`--max-reactions-per-cycle`、`--reaction-cutoff`（既定4 Å）、`--max-reaction-cutoff`（既定6 Å）、
`--fixed-reaction-cutoff`、`--reaction-minimization-iterations`、`--reaction-minimization-radius`（既定12 Å、0で系全体）、
`--post-reaction-md-steps`、`--relaxation-stages`（JSON文字列かJSON/YAMLファイル。書式は
[crosslinking/README.md](../crosslinking/README.md#反応距離の上限と反応後の緩和)）、`--stop-after-idle-cycles`（既定200）、
`--template-workers`（既定4）、`--final-relaxation-steps`（既定1000）、`--final-md-steps`（既定1 ns）、
`--save-trajectory`。`--remove-sol-fraction`で、現像を模して周期次元3の網目だけを残し、
未反応モノマーと網目につながっていないオリゴマー（ゾル分）を除いた構造を物性評価に使います
（ゾル分率・削除数・ゲルだけで数えた環の反応率`gel_epoxide_conversion`を記録。ゲル化点未満ならそのレプリカは停止）。
溶媒は`--solvent NAME=SMILES`（複数可）と`--solvent-molecules`で架橋プロジェクトに指定します。
YAMLや段の設定は作成時に`project.json`へ取り込まれ、初回run後に設定が変わると停止します。
条件を変える場合は新しいプロジェクトを作成してください。

`--activation-count 0`（`n0`）は開始点なし・反応なしで、単分子の★構造をそのまま熱履歴と300 K物性に回します。

### 保存済みネットワーク（`network`）

```bash
python $NP network workflows/network_properties/output/saved --network final.structure.json --solvent water=O
```

入力JSONをproject内にコピーし、各レプリカは熱履歴から開始します（共通トポロジー、熱履歴だけが独立）。

## ローカルGPUへレプリカをジョブ投入する

`--replicates N`で指定した本数を、`submit`でこのマシンのGPUへ投入できます。
GPU番号は`nvidia-smi`に表示される物理GPU番号です。
単分子プロジェクトを指定すると、その`crosslinks/`の架橋プロジェクトもまとめて投入します。
架橋レプリカNは単分子レプリカNの完了を待って開始し、待っている間は他のジョブにGPUを使います。
単分子レプリカNが失敗すると、架橋レプリカNは`blocked`になります（`reseed`で単分子を直してから再度`submit`）。
`--monomer-only`を付けると、`crosslinks/`があっても単分子のレプリカだけを投入します
（単分子を全部終わらせてから架橋を流す場合。架橋は後で同じディレクトリを`submit`すれば、完了済みの単分子を省いて投入されます）。
`--mechanical`を付けると、300 K構造ができている架橋・ネットワークのレプリカの弾性計算だけを投入します。

```bash
# 例：4レプリカ、base seed 31415の単分子と、2つの架橋密度
python workflows/network_properties/network_properties.py monomer \
  workflows/network_properties/output/epoxy_gpu --replicates 4 --seed 31415 --platform CUDA
python workflows/network_properties/network_properties.py crosslink \
  workflows/network_properties/output/epoxy_gpu --activation-fraction 0.05 0.1 --max-reactions-per-cycle 50

# GPU割り当てとseedを確認する（このコマンドでは起動しない）
python workflows/network_properties/network_properties.py submit \
  workflows/network_properties/output/epoxy_gpu --gpus 0,2 --dry-run

# GPU 0とGPU 2に投入。各GPUで1ジョブずつ実行し、残りは空き次第実行
python workflows/network_properties/network_properties.py submit \
  workflows/network_properties/output/epoxy_gpu --gpus 0,2

python workflows/network_properties/network_properties.py status \
  workflows/network_properties/output/epoxy_gpu
```

GPUが4台あれば`--gpus 0,1,2,3`で4本を並列実行できます。
GPUが1台なら`--gpus 0`で全レプリカを順番に実行します。
番号を省略すると、このプロセスから利用可能な全GPUを使います。
既存の`CUDA_VISIBLE_DEVICES`で制限されたGPUの範囲を超えて割り当てません。
内部ではGPU UUIDを使って各ジョブから1台だけを見せるため、物理GPU番号の取り違えを防ぎます。

各レプリカのseedは`project.json`の`seed + replica番号（1始まり）`です（`reseed`した場合を除く）。
上の例では31416, 31417, 31418, 31419になります。
単分子ではpack・熱履歴・Tの物性へ、架橋では開始点の選択・架橋反応・熱履歴・物性評価へこのseedを渡します。
単分子レプリカ1本（構築からTの物性まで）、架橋レプリカ1本（架橋から溶媒評価まで）がそれぞれ1ジョブです。
`network`プロジェクトでは構造が共通で、熱履歴だけが独立になります。

`submit`はバックグラウンドの管理プロセスを起動して戻ります。
端末を閉じてもジョブは継続し、`jobs/replica_01.log`などに標準出力・エラーを保存します。
管理プロセスのログは`jobs/launcher.log`、PIDは`jobs/launcher.json`、
GPU・seed・ジョブ状態は`jobs/plan.json`と`jobs/status.json`に記録します。
各レプリカにも`execution.json`を保存します。

### seedを変えて失敗したレプリカだけ再計算する

seedによってはMDが発散するなどしてレプリカが失敗します。失敗したレプリカだけに新しいseedを
割り当て、そのレプリカだけを最初から計算し直せます。

```bash
# 新しいseedを確認（何も変更しない）
python workflows/network_properties/network_properties.py reseed \
  workflows/network_properties/output/epoxy_gpu --replica 3 --dry-run
# 失敗したreplica 3とreplica 4に新しいseedを割り当て、そのままGPU 0へ投入
python workflows/network_properties/network_properties.py reseed \
  workflows/network_properties/output/epoxy_gpu --replica 3 --replica 4 --submit --gpus 0
```

新しいseedは`seed + k×レプリカ数 + replica番号`（k = 1, 2, …）のうち、そのprojectでまだ
使っていない最小の値です（`--seed`で1レプリカに明示も可）。初回のseedや他のレプリカと重なりません。
旧seedの出力は`failed_attempts/replica_03_seed31418_1/`へ移し（エラーログも含めて保存）、
レプリカは最初（単分子ならpacking、架橋なら単分子構造への開始点作成）からやり直します。
単分子レプリカをreseedすると、そのレプリカから作った架橋レプリカは出発構造が変わるため実行を拒否します
（架橋側も同じ番号をreseedしてください）。各段階がseedに依存するため、途中からは再開しません。
割り当ては`seed_overrides.json`に記録し、`project.json`は変更しません
（他のレプリカの実行条件の指紋はそのままです）。`status`、`report.md`、`replica_seeds`にも新しいseedが出ます。
対象は`failed`状態のレプリカだけです。完了済みや、強制終了で状態が残っていないレプリカは`--force`が必要です。
実行中のレプリカやGPUジョブの投入中は拒否します。`--submit`を付けない場合は、あとで`submit`すれば
未完了のレプリカとして実行されます。

`status`は単分子と各架橋プロジェクトのレプリカごとに、現在の計算段階（単分子の構築、架橋、熱履歴の各温度、
Tの段階、NVT、弾性、溶媒ごとのpack・挿入自由エネルギー・拡散MD）と経過時間、架橋中なら最新のcycleと反応率、
ジョブログの最終行を表示します（投入前でも段階の記録から表示します）。段階は`replica_XX/progress.json`に記録し、直近の履歴も残します。
ランチャーやワーカーのPIDが存在しない、または終了済みで回収されていない（zombie）場合は
`status`が警告します。その場合は各レプリカの状態を確認し、再度`submit`すると未完了分を再開します。

計算がエラーで終了すると、失敗した段階とPythonのtracebackを`replica_XX/error.log`へ書き出します。
架橋構築で失敗した場合は、`crosslinking.py`のエラー行と標準エラー出力の末尾も含めます。
ランチャーは原因の1行を`jobs/status.json`の`error`へ記録し、今回の試行分のログ末尾300行を
`jobs/replica_XX.error.log`へ保存します。SIGKILLなどのシグナルで強制終了された場合は
シグナル名を記録します（SIGKILLはメモリ不足時のOOM killerでも発生します）。
単分子の構築結果は`replica_XX/build/results/monomer/`、架橋の結果は`replica_XX/crosslinking/results/network/`、
架橋の出発構造・到達反応率・構造ハッシュは`replica_XX/network_origin.json`に記録します。

同じprojectへの二重投入、および同じレプリカの二重実行はロックで拒否します。
失敗したレプリカがあっても他のレプリカは継続し、全ジョブ終了時にレポートとCSVを集計します。
`completed`は計算処理の正常終了を示し、物性値の収束判定とは別です。
再度`submit`すると完了済みレプリカを省略し、未完了レプリカを再開します。
完了済みの架橋構築はハッシュを検査して再利用します。
架橋構築そのものには途中再開機能がないため、中断した構築は同じseedで最初からやり直します。
その際は当該レプリカの未完了構築出力を作り直し、局所Sageテンプレートのキャッシュは再利用します。
最大cycleまでに反応率の目標へ届かなくても、その実現構造で後続評価を行います。
到達反応率と未達フラグをレポートに示し、未達レプリカはチェック通過値の集計に含めません。
投入の単位（単分子とその架橋、または1つのプロジェクト）が異なれば、同じGPUを指定しても共有は制限しないため、
割り当ては投入間で調整してください。

単一レプリカを手動で実行することもできます。

```bash
CUDA_VISIBLE_DEVICES=0 python workflows/network_properties/network_properties.py run \
  workflows/network_properties/output/epoxy_gpu --replica 1
```

通常の`run`は全レプリカを順次実行します。単分子プロジェクトに`--all`を付けると、続けて全架橋プロジェクトも実行します。
集計ファイルにもロックをかけ、並列ジョブや手動`report`との書き込み競合を防ぎます。
`--smoke --replicates 2`のように、動作確認でも複数レプリカを明示できます。

## 熱履歴と300 K物性

| 段階 | 既定条件・出力 |
|---|---|
| 初期緩和と加熱 | 最小化、300 Kから100 K刻みで加熱、各100 ps |
| 高温緩和 | 800 K、NPT、1 ns |
| 冷却 | **800→100 K**、20 K刻み、各500 ps、平均40 K/ns |
| 冷却解析 | 各温度plateauの後半を使用。密度の双曲線フィットからTg、線・体積CTE（フィット曲線のTg∓50 Kでの傾き） |
| 温度復帰 | 100→300 K、20 K刻み、各100 ps |
| 300 K再緩和 | NPT、1 ns |
| 乾燥NPT観測 | 300 K、1 atm、10 ns、1 ps間隔 |
| 乾燥NVT観測 | 冷却・復帰後の固定セルで1 ns再緩和、10 ns観測 |
| 独立seed | 既定では架橋構築から独立した4レプリカ。保存構造入力時は共通トポロジー |

熱履歴の終了時に`replica_XX/tg_density_volume.png`を書きます。冷却中の温度–密度・温度–体積（サンプルと10 K平均）、Tgフィット曲線と、CTEを取るTg∓50 Kでの接線、Tg、昇温各段階の終了時の密度・体積を1枚に描きます。昇温中の点が冷却曲線から大きく外れていれば、高温での沸騰・空洞化を疑ってください（`report`でも既存の結果から作り直し、report.mdに載せます）。

冷却は温度plateauを使う段階冷却です。800 Kにも観測用plateauを置くため、
冷却区間の合計は既定18 nsです（800→100 Kの温度差に対応する分は17.5 ns）。
時間は刻み幅とタイムステップに丸め、各段階の`done.json`に実時間を記録します。
既定タイムステップは1 fs（X–H結合を拘束。熱履歴・300 KのNPT/NVT・溶媒で共通）。
2026-11-12以前の既定は0.5 fsで、時間（ps）は同じまま刻みを倍にしたので、step数は半分です。
体積の試行間隔も50 step（50 fs、以前は0.5 fs×100 step）にして同じ頻度に保っています。
弾性だけはX–H拘束を外し、別に0.25 fsを使います（`timestep_fs`の影響を受けません）。
体積制御はwheelのconstraint group座標に対する
周期架橋対応Monte Carlo更新です。一般的な分子単位barostatへの置き換えはしません。

**100 Kまで冷却した後に300 Kへ戻した構造**を、弾性・溶媒評価に使います（単分子ではTの段階もここから続けます）。
冷却途中の300 K構造も別名で保存します。再加熱でTgを再評価する設定ではありません。
800 Kで化学的に分解する場合でも、固定トポロジーMDには分解機構は含まれません。

**未架橋・低反応率の系では高温側を下げてください。** 反応率0のTMP triglycidyl ether（3000原子）では、
800 K・1 atmで密度が300 Kの0.50倍（1.07→0.54 g/cm³）まで下がり、セル内に半径6 Åを超える空洞が
できました（沸騰・キャビテーション）。この空隙は冷却後のガラスに残り、Tg/CTEのフィットと300 Kの揺らぎを
壊します。網目は蒸発できないので、これは単分子・オリゴマーに特有です。
単分子プロジェクトも既定は800 Kのままです。`--high-temperature 600`などで最高温度を変えられます（`project.json`の`high_K`）。
各段階の終了時密度は`dry_result.json`の`thermal_history_density`に記録し、高温緩和後の密度が
開始時の0.7倍を下回ると`possible_vaporization: true`として、冷却フィットと300 K NPTのチェックを不合格にします。

RadonPyと同じ熱力学的揺らぎ関係を使う出力は次の通りです。

| 物性 | 関係・単位 |
|---|---|
| 密度 | 質量 / 体積、g/cm³ |
| 定圧比熱Cp | Var(H) / (kB T² m)、J/(kg K) |
| 等温圧縮率κT | Var(V) / (kB T ⟨V⟩)、1/Pa |
| 等温bulk modulus KT | 1/κT、GPa |
| 体積CTE αV | Cov(V, U+PV) / (kB T² ⟨V⟩)、1/K（= Cov(V,H)/…、Kは体積と無相関） |
| 線CTE αL | αV / 3（等方性を仮定）、1/K |
| 定積比熱Cv | Cp − T αV² / (ρ κT)、J/(kg K) |
| 断熱圧縮率κS | κT Cv / Cp、1/Pa |
| 断熱bulk modulus KS | 1/κS、GPa |
| bulk音速 | √(KS/ρ)、m/s。固体の縦波音速とは異なる |
| 直接Cv | NVTでのVar(E) / (kB T² m)、J/(kg K) |
| 固定電荷モデルの誘電応答 | NVTの連続的電荷変位から誘電テンソル、等方平均 |
| 幾何補助量 | wheelのpoint-probe空隙率など。溶媒アクセス可能体積とは異なる |
| ネットワーク形状 | 連結成分、最大成分の質量割合、cycle rank、周期セルを貫く結合cycleの有無 |

Hは`U + K + Pext V`です。体積CTEの共分散には配置エンタルピー`U + Pext V`を使います。NPTでは運動量が配置・体積と独立でCov(V, K) = 0のため期待値は同じで、Kを含めると雑音だけが増えます（2026-11-12以前はKを含めていたため、短い観測でCTEが負になりやすかった）。比熱は古典・拘束付きSageモデルの値で、量子補正は行いません。
5ブロックの誤差、体積・エネルギーの有効サンプル数と前後半drift、
replica間の標本SDとSEMを出力します。1 replicaだけの場合、replica間誤差はnullです。
ガラス・ゴムの線CTEは、双曲線フィットの漸近線ではなく、フィット曲線そのもののTg−50 K／Tg+50 Kでの傾きから
αL = −(dρ/dT)/(3ρ(T)) で求めます（`cte_offset_K`、既定50 K）。漸近線は遷移幅が広いとデータ範囲外への外挿になり、
曲線がどこでも右下がりなのにガラス側が負のCTEになることがありました。Tg±50 Kが冷却範囲の外ならフィットのチェックは不合格です。
遷移幅wが広いとガラス側にゴム側の傾きが混ざります（混ざる割合`cte_mixing_fraction` = (1 − 50/√(50²+w²))/2）。
既存の結果は`report`／`collect`の時に`cooling.csv`から再フィットし、`dry_result.json`とTg図を更新します（MDの再実行は不要）。
Tgは冷却速度依存の見かけ値です。フィットのR²、遷移幅、端への張り付き、
単純な直線よりAICが十分改善するかも確認します。

300 Kの点プローブ空隙率はNPT観測終了後のスナップショットに対して計算し、
`dry_result.json`の`structure.point_probe_void_fraction`とレプリカ比較CSVに保存します。
本計算のランダム点数は10000点です。0〜1の割合で、時間平均・有限サイズ溶媒の
アクセス可能体積・空隙ネットワークの連結性を表すものではありません。

誘電応答はPMEの導体境界条件を仮定した固定電荷モデルの応答です。
周期境界のbranch jumpを避けるため、各原子の最小像変位を累積して
全電荷変位を作ります。全体は中性である必要があります。
電子分極は含まれず、ガラスの配向緩和が観測時間より長ければ収束しません。
そのため、誘電物性を自動的に「収束済み」とは判定しません。

## 網目のグラフ解析

300 KのNPT観測終了後の構造について、化学グラフから網目の構造量を求め、
`dry_result.json`の`network`とレプリカ比較CSVに保存します（MDは不要で、2万原子で1秒以下）。
保存済みの構造だけを解析する場合は次を実行します。

```bash
python workflows/network_properties/network_graph.py replica_01/relaxed_300K.structure.json \
  --histogram-csv strands.csv --histogram-png strands.png   # 1鎖1行のCSVと分布の図
```

網目鎖は、もとのモノマーの区切りに関係なく、架橋点から次の架橋点までを1本と数えます
（重合でつながった複数モノマーにまたがる鎖も1本）。鎖ごとに次を求め、全分布を出力します。
- 長さ`strand_heavy_atoms_histogram`：架橋点の間にある重原子数。途中の環系は側の原子も含めて全部数えます
  （フルオレンの9位を通る鎖ではフルオレン13原子をすべて含む）。
- 分子量`strand_molecular_weight_*`：同じ原子にHを含めた分子量（g/mol）。ヒストグラムは25 g/mol幅。
- 主鎖の結合数`strand_backbone_bonds_histogram`：架橋点の原子から架橋点の原子までの最短経路の結合数（参考）。

`strand_table`には鎖1本ごとの3つの値が入ります。reportでは各レプリカの`network_strands.csv`（1鎖1行）と
`network_strands.png`（3つの分布）にも出力します。鎖の平均長・平均分子量、鎖と架橋点の本数は
レプリカ比較とディレクトリ横断集計にも入ります。

H原子を除いた重原子グラフで、縮環を含む3〜8員環の環系を1つの節点にまとめます
（鎖の途中の芳香環は通過点で、分岐点とは数えません。周期境界を通って閉じる環は環とみなしません）。

| 出力 | 定義 |
|---|---|
| `sol_mass_fraction` | 周期境界をまたいで自分自身の像とつながらない分子（未反応モノマー・オリゴマー）の質量分率 |
| `dangling_mass_fraction` | ゲルのうち、1本の結合（橋）だけでつながり、その先に無限につながった部分がない部分の質量分率 |
| `elastically_ineffective_loops`、`loop_mass_fraction` | 同じ架橋点に戻るループの数と質量分率。ループを除いて架橋点を数え直し、これを繰り返すので、付け根の節点は結合が2本に減れば鎖の一部に戻ります |
| `active_mass_fraction` | ゾル、ダングリング、ループを除いた有効網目の質量分率 |
| `junctions`、`junction_functionality` | ループ除去後の有効網目で、有効な結合を3本以上持つ節点（架橋点）とその官能度の分布 |
| `effective_strands`、`strand_density_nu_e_*` | 架橋点どうしを結ぶ鎖の数と、その数密度ν_e |
| `Mc_g_mol` | ρN_A/ν_e。実質は「セル内の全質量÷有効鎖本数」の見かけの値で、ゾルなどの質量も含みます |
| `Mc_active_network_g_mol` | 有効網目の質量÷有効鎖本数。架橋点の原子を含みます |
| `strand_molecular_weight_g_mol` | 鎖そのものの分子量の統計。架橋点の原子を含みません |
| `cycle_rank_density_mmol_cm3` | ξ/V（ξ = 鎖の数 − 架橋点の数） |
| `rubbery_shear_modulus_estimates_MPa` | ゴム状態の弾性率の目安。affine：ν_e kT、phantom：ξkT/V（300 K。からみ合いは含まない） |
| `polar_groups_*` | 脂肪族／フェノール性OH、エーテルO、未反応エポキシ環・オキセタン環、カルボニルO、N（H数別）、芳香環、水素結合の供与H・受容原子、反応活性点（`cationic_activity`）の数・mmol/cm³・mmol/g |

水素結合の受容原子はNとOを官能基で判定します（RDKitのLipinski定義に準拠。アミド型N、
酸のOH、ニトロのN、ピロール型の芳香族Nなどは数えません）。供与は原子ではなくHの本数を数えます。
OH*の代理表現では、活性な鎖末端もOHとして数えます（`active_cationic_sites`と同数）。
結果には`method_version`が入り、定義が古い値は比較表に残りますがチェック通過には数えません。
架橋エポキシはガラス状態なので、ν_eやゴム状態の弾性率はFlory–Rehner式など膨潤理論の入力用の値で、
300 Kのガラスの弾性率（「弾性」の節）とは別物です。

## OVITOでの可視化

各段階の`state.xml`はOpenMMの状態（座標・速度・セル、nm）だけで、元素・結合・電荷を含みません。
それらは同じディレクトリの`final.structure.json`（PAVES形式、Å）にあります。
`to_ovito.py`はこのJSONから、OVITOでそのまま開けるLAMMPS data（atom_style full：元素名の原子タイプ、
部分電荷、分子番号、結合、セル内へwrapした座標とimage flag）を書きます。

```bash
# 1構造 -> relaxed_300K.data
python workflows/network_properties/to_ovito.py replica_01/relaxed_300K.structure.json
# 元素・結合はJSON、座標とセルはOpenMMのstate.xmlから
python workflows/network_properties/to_ovito.py replica_01/02_high_equilibration/final.structure.json \
  --state replica_01/02_high_equilibration/state.xml -o high.data
# レプリカの全段階：段階ごとの.dataと、段階順の多フレームextended XYZ（セルの変化を含む）
python workflows/network_properties/to_ovito.py replica_01 --trajectory replica_01/history.xyz
```

extended XYZには結合が入りません（OVITOのCreate bondsで付けます）。各フレームの段階名は`phase`属性に入ります。
numpyだけで動くので、OVITOの環境でも実行できます。

## 弾性

弾性は既定では計算しません。架橋膜（または`network`）の300 K構造ができた後に、
`mechanical <プロジェクト> [--replica N]`か`submit <プロジェクト> --mechanical`で実行します。
結果は`replica_XX/mechanical/`に保存し、`report.md`の「Elasticity」と架橋シリーズのYoung率に出ます。
単分子の液体には実行できません。

弾性は300 K構造から6成分×正負1%のひずみへ分岐します。
各分岐で200 ps NVT緩和、1 ns観測を行い、全ポテンシャルの微小セル変形差分
（PMEも含む）と運動項から応力を計算します。Voigt順序はxx, yy, zz, yz, xz, xy、
shearは工学ひずみです。
2万原子の網目では1 psごとの応力の揺らぎが30〜70 MPaあり、0.5%のひずみ・短い観測では
応答（約25 MPa）が揺らぎに埋もれて行列が非対称・不安定になりました。既定の1%・1 nsでは
弾性定数の標準誤差は伸縮成分で約0.1 GPa、ずり成分で約0.05 GPaの見込みです
（1 psごとの応力をほぼ独立とみなした概算）。12分岐×1.2 ns（0.25 fs）で、
RTX 5070 Ti 1枚で1レプリカ約1日かかります。

出力は6×6応力–ひずみ接線行列、ブロック誤差、対称性のずれ、
Voigt/Reuss/Hillのbulk/shear modulus、HillのYoung率・Poisson比です。
行列に非正の固有値があれば等方弾性率を確定しません。
仮想ひずみに拘束のvirialを欠落させないため、弾性計算ではX–H結合を拘束せず、
Sageの調和結合項として読み込み（`load_network(..., constrained=False)`）、自然質量・0.25 fsを使います。
NPTの拘束付きKTとの違いにはこのモデル差も含みます。
以前は拘束付きの系から拘束だけを外していたため、X–H結合に結合項がなく（OpenFFは拘束した結合に
結合項を書きません）、ひずみ0でも約−0.7〜−0.8 GPaの見かけの応力が出ていました。
2026-10-30以前の弾性結果は再利用も集計もされません。`mechanical/provenance.json`の手法バージョン、
入力構造のハッシュ、ひずみと時間の設定が一致しない結果は`mechanical.stale-N/`へ移して再計算します。
wheelの更新は自動では検出しません。応力が変わる変更では`mechanical.py`の`METHOD_VERSION`を上げてください。

弾性は、ひずみ幅、応力の標本誤差、セル長依存の
比較が必要なため、自動収束フラグを常にfalseにしています。
値は比較検証用の推定値として保存します。

## 溶媒和自由エネルギーと拡散係数

各溶媒について、同じ`relaxed_300K.structure.json`から2つに分岐します。

1. **乾燥樹脂への希薄な試験分子挿入**：N=0の樹脂へ、最初は非相互作用の1分子を置き、
   soft-core LJと電荷を段階的にcoupleします。既定22窓、各100 ps緩和＋500 ps観測。
   wheelのMBAR解析から`coupling_free_energy_kJ_mol`を得ます。
   これは対応する孤立した柔軟分子から樹脂環境へ移す過剰自由エネルギーです。
   `chemical_potential_relative_kJ_mol`は有限セルの濃度項を含む別の量として保存します。
   参照数密度は0.602214076 nm⁻³（1 mol/L）で、`standard_density_nm3`で変更できます。
2. **有限溶媒量での自己拡散**：指定数をpack、最小化、100 ps NVT、1 ns NPT緩和、
   10 ns NPT観測。既存wheelのCOM変位追跡で周期境界、樹脂並進、セルのaffine変形を除き、
   複数時間原点のMSDと`D = slope(MSD)/6`から評価します。

自由エネルギーには窓間overlap、有効サンプル数、前後半の差、誤差、lambda経路を保存します。
1つのcoupling方向とseedだけで精度を保証せず、独立初期位置・逆方向・窓追加・セルサイズを
比較してください。特に高密度ガラスでは異なるcavity間の移動が遅い場合があります。
液体溶媒から樹脂への移行自由エネルギーが必要なら、同じモデルによる液体側の
過剰化学ポテンシャルも求め、その差を取る必要があります。ここでは気相／孤立分子基準です。

拡散解析では最低観測時間、MSD変位量、log(MSD)–log(t)の傾きなどを確認します。
拡散域に届かない場合は`diffusion_coefficient_m2_s`をnullにし、
候補の傾きだけを別フィールドに保存します。有限サイズ補正は行いません。
自由エネルギーは希薄挿入、拡散は指定した有限濃度です。両方を希薄条件に近づけるには
`--solvent-molecules 1`を使用し、セルサイズ依存も確認してください。

指定した溶媒数は**平衡吸収量ではありません**。十分な空間へpackできなければ停止し、
セルを勝手に広げません。膨潤は溶媒MDのNPTで生じます。

## 出力と再開

```text
project.json / run_signature.json
seed_overrides.json                 # reseedしたレプリカの新しいseedと履歴
failed_attempts/                    # reseed前の失敗した試行（レプリカとしては集計しない）
definitions/                        # 単分子・架橋：コピーした開始点・反応ルール
input.structure.json                # networkのみ：共通入力構造
summary.json / report.md
replica_properties.csv              # 各物性を行、各replicaの値・確認状態を列にした一覧
property_summary.csv                # 平均・標本SD・SEM・有効件数・チェック通過件数
crosslink_series.csv / .png         # 単分子のみ：架橋密度ごとの比較
jobs/                               # submit時：投入計画・状態・レプリカごとの実行ログ
replica_01/
  execution.json                    # レプリカ番号・seed・実行時に見えるGPU
  build/results/monomer/            # 単分子：pack・緩和したモノマー構造（build/done.jsonに記録）
  crosslinking/results/network/     # 架橋：このレプリカ専用の架橋構造・履歴・設定
  network_origin.json               # 架橋：出発構造（単分子レプリカ）・反応率・原子数・構造ハッシュ
  cooling.csv
  cooling_300K.structure.json
  relaxed_100K.structure.json
  relaxed_300K.structure.json
  dry_result.json
  00_... / ... / 06_production_300K/  # samples, state, structure, done marker
  07_nvt_300K/                      # direct Cv and fixed-charge dielectric response
  363.15K/                          # 単分子：Tの段階（result.json、relaxed_363.15K.structure.json）
  mechanical/result.json            # 任意
  solvents/water/
    packed/
    infinite_dilution_free_energy/  # windows, overlap, checkpoints, result
    diffusion/                     # MSD, diagnostics, OVITO frames, final structure
    result.json
crosslinks/af0.05/ …                # 単分子のみ：架橋プロジェクト（上と同じ構成）
```

同じ`run`コマンドで再開します。熱履歴は完了済みの温度・段階のstateから戻し、
中断した段階を初めから実行します。弾性は手法と入力が同じ場合に限り、完了した正負分岐を
再利用します。launcherが強制終了するとワーカーも止まり、まだ動いているレプリカがあれば`submit`は拒否します。
`report`と`collect`は、実行開始後に`project.json`を変更したプロジェクトを集計しません。
溶媒MD・自由エネルギーはwheelのcheckpointから再開します。
熱履歴の再開はOpenMMの乱数状態まで一致する厳密再現ではありません。
`network`で共通構造を入力した場合は、異なる架橋トポロジーによるばらつきを含みません。

`report.md`では各物性をreplicaごとに横並びにし、平均・標本SD・SEM、値が得られた本数、
既存の収束・フィットチェックを通過した本数も表示します。
乾燥熱物性、Tg/CTE、弾性、溶媒ごとの挿入自由エネルギー・拡散係数を
同じ形式で集約します。誘電テンソル・弾性テンソルの各成分もCSVに含めます。

`replica_properties.csv`のstatus列は`checks_passed`、`unchecked`、`smoke`、
`pending`（未完了）、`unavailable`（結果はあるが値が得られない）のいずれかです。
未計算・null・非有限値は空欄とし、平均の計算に0として加えません。
集計CSVの`mean/std/sem`は有限な候補値すべて、`checked_mean/checked_std/checked_sem`は
チェック通過値だけの集計です。レポートでは未確認値と、それを含む平均に`*`を付けます。
SD/SEMはreplica間のばらつきで、1本しか値がない場合は未定義（空欄／null）です。
各計算内のMBAR誤差やblock SEMは、元のreplica別JSONに保存しています。

集計は`run`で各replicaが完了するたびに更新します。既存の計算結果にも
`report`コマンドで適用でき、MDの再実行は不要です。
途中の結果を手動集計する場合、溶媒の自由エネルギー・拡散の個別`result.json`も読み込みます。

屈折率、Abbe数、電子応答にはQM分極率などが必要です。
また、周期的に連結した架橋の凝集エネルギー密度・Hildebrandパラメータには
一意な孤立分子基準がありません。これらは未計算とし、summary.jsonに理由を記録します。

## 動作確認

```bash
export OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
NP=workflows/network_properties/network_properties.py
python $NP monomer workflows/network_properties/output/smoke --platform CPU --smoke --replicates 2
python $NP crosslink workflows/network_properties/output/smoke --activation-fraction 0.1 0.2 \
  --solvent water=O --solvent-molecules 1
python $NP run workflows/network_properties/output/smoke --all
python $NP mechanical workflows/network_properties/output/smoke/crosslinks/af0.1

python -m unittest discover -s tests -v
PAVES_NETWORK_INTEGRATION=1 python -m unittest discover -s tests -p test_network_properties.py -v
```

**`--smoke`は全段階を極端に短縮する接続テストです。値に物理的意味はありません。**
reportにも明記し、収束済みにはしません。

## 参照

- [RadonPy公式：対象物性と評価手順](https://github.com/RadonPy/RadonPy)
- [Hayashi et al., RadonPy, npj Comput. Mater. 8, 222 (2022)](https://doi.org/10.1038/s41524-022-00906-4)

RadonPyを依存として実行するのではなく、物性範囲と揺らぎ式を参照し、
PAVESが保存したSageネットワーク上で評価します。
