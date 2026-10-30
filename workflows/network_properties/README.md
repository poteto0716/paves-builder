# 架橋ネットワークの熱・力学・溶媒物性評価

既定では、各レプリカがモノマーのpack・開始点作成・架橋反応から独立に実行し、
自分で生成した構造を熱緩和・乾燥物性・溶媒評価へ渡します。
架橋後の結合、部分電荷、安定原子IDを保持し、その後は固定トポロジーのSage/OpenMMで計算します。
`--network final.structure.json`を明示すると、保存済みの共通ネットワークから評価できます。
既存のPCFF鎖構築用`thermomechanical.py`とは別の入力経路です。
周期境界をまたぐ架橋を有限分子MOL2として再importしません。

## 実行

架橋構築と同じ[Python 3.11 / OpenFF環境](../crosslinking/environment.yml)と
公開wheel `dist/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl`を使用します。
リポジトリのルートで実行してください。

```bash
conda activate paves-reactive
python workflows/network_properties/network_properties.py new \
  workflows/network_properties/output/epoxy \
  --replicates 4 --seed 31415 --max-reactions-per-cycle 50 \
  --solvent water=O --solvent methanol=CO --solvent 'acetone=CC(=O)C' \
  --solvent-molecules 10 --platform CUDA

# 必要に応じて、作成されたproject.jsonで計算時間やseedを変更する。
python workflows/network_properties/network_properties.py run \
  workflows/network_properties/output/epoxy
python workflows/network_properties/network_properties.py report \
  workflows/network_properties/output/epoxy
```

溶媒を省略すれば乾燥物性だけを実行します。溶媒は複数指定可能です。
Sageが割り当て可能な単一の中性分子が対象で、イオン、ラジカル、混合SMILESは
既存の溶媒APIで拒否されます。水もSageのモデルです。
溶媒数は各溶媒で共通の既定値を使います。異なる濃度の比較は別projectで作成します。
SMILESが同じでも温度、溶媒数、樹脂の架橋率が違えば、別の条件の物性です。

`--skip-mechanical`で追加MDの多い弾性評価を無効にできます。
熱伝導率（RNEMD）は架橋ネットワークでMDが発散したため、このワークフローから外しました。
作成済みのproject.jsonにもすべての設定を記録します。
架橋構築の既定は総原子数20000を超える最小分子数、目標反応率99.5%、最大5000 cycle、
トラジェクトリ保存OFFです。モノマーはTMP triglycidyl etherです。
`--remove-unreacted-monomers`で、架橋後に他の分子と結合しなかったモノマーを
除いた構造を物性評価に使います（既定OFF。抜けた体積は熱履歴の昇温NPTで緩和）。
packing（既定0.5 g/cm³、`--packing-density`。PAVESの局所緩和あり、環を貫通する結合があれば停止）後、1 atmの初期NPTの前に300 K・200 MPa・50 psの加圧NPTでvoidを潰します。
`--compression-pressure-mpa`、`--compression-steps`（1 fs刻み、0で無効）で変更できます。
`--epoxy-smiles`、`--total-atoms`、`--target-conversion`、`--cycles`、
`--activation-fraction`または`--activation-count`で変更できます。
反応距離の上限は`--max-reaction-cutoff`（既定7 Å）で変えられます。
開始距離は`--reaction-cutoff`、距離の固定は`--fixed-reaction-cutoff`です。
反応後の緩和は`--reaction-minimization-iterations`、`--reaction-minimization-radius`（反応点から何Å以内を最小化するか。既定12、0で系全体）、`--post-reaction-md-steps`、
`--relaxation-stages`（JSON文字列かJSON/YAMLファイル）で設定します。
段の書式は[crosslinking/README.md](../crosslinking/README.md#反応距離の上限と反応後の緩和)を参照してください。
例:

```bash
python workflows/network_properties/network_properties.py new workflows/network_properties/output/epoxy_all \
  --replicates 4 --seed 31415 --max-reactions-per-cycle 50 --solvent water=O --platform CUDA \
  --max-reaction-cutoff 5 --post-reaction-md-steps 1000 --reaction-minimization-iterations 5000 \
  --relaxation-stages workflows/crosslinking/relaxation.example.yaml
```

反応終了後は0.1 fs×1000 stepの緩和と1 atm・1 nsのMD（`--final-relaxation-steps`、`--final-md-steps`）を経て最終構造を保存します。
YAMLファイルの中身は`new`の時点で`project.json`に取り込まれます。そのため、後でファイルを書き換えても既存projectには影響しません。

活性点が網目に埋もれて反応が止まった後の空回りは、`--stop-after-idle-cycles`（既定200、0で無効）で打ち切ります。
未知の反応テンプレートは`--template-workers`（既定4）個のプロセスで並列に計算します。
詳細は[crosslinking/README.md](../crosslinking/README.md#計算速度)を参照してください。
これらは`project.json`の`crosslinking`（`max_reaction_cutoff`、`minimization_iterations`、
`post_reaction_md_steps`、`relaxation_stages`）にも記録されます。初回run前なら、そこを直接編集しても構いません。
独自の反応ルールは`--activation-yaml`、`--reactions-yaml`で渡します。
これらはproject内の`crosslinking/`へコピーし、SHA-256で変更を検査します。
`project.json`の`crosslinking`以下には全構築条件を保存します。
初回run後に設定が変わると既存結果との混在を防ぐため停止します。
条件を変える場合は新しいprojectを作成してください。

保存済み構造を共通入力として使う場合だけ、作成時に`--network`を指定します。
その場合は入力JSONをproject内にコピーし、各レプリカは熱緩和から開始します。
`--network`とエポキシ構築条件は同時に指定できません。

## ローカルGPUへレプリカをジョブ投入する

`new --replicates N`で指定した本数を、`submit`でこのマシンのGPUへ投入できます。
GPU番号は`nvidia-smi`に表示される物理GPU番号です。

```bash
# 例：4レプリカ、base seed 31415のプロジェクトを作成
python workflows/network_properties/network_properties.py new \
  workflows/network_properties/output/epoxy_gpu \
  --replicates 4 --seed 31415 --max-reactions-per-cycle 50 \
  --solvent water=O --platform CUDA

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

各レプリカのseedは`project.json`の`seed + replica番号（1始まり）`です。
上の例では31416, 31417, 31418, 31419になります。
既定では、モノマーpack・開始点の選択・架橋反応・熱緩和・物性評価へこのseedを渡します。
各レプリカの架橋ネットワークを別々に生成し、構築から溶媒評価までの一連を1ジョブとして実行します。
`--network`モードでは架橋構造が共通で、熱履歴だけが独立になります。

`submit`はバックグラウンドの管理プロセスを起動して戻ります。
端末を閉じてもジョブは継続し、`jobs/replica_01.log`などに標準出力・エラーを保存します。
管理プロセスのログは`jobs/launcher.log`、PIDは`jobs/launcher.json`、
GPU・seed・ジョブ状態は`jobs/plan.json`と`jobs/status.json`に記録します。
各レプリカにも`execution.json`を保存します。

`status`は各レプリカの現在の計算段階（架橋構築、熱履歴の各温度、NVT、弾性、
溶媒ごとのpack・挿入自由エネルギー・拡散MD）と経過時間、架橋中なら最新のcycleと反応率、
ジョブログの最終行を表示します。段階は`replica_XX/progress.json`に記録し、直近の履歴も残します。
ランチャーやワーカーのPIDが存在しない、または終了済みで回収されていない（zombie）場合は
`status`が警告します。その場合は各レプリカの状態を確認し、再度`submit`すると未完了分を再開します。

計算がエラーで終了すると、失敗した段階とPythonのtracebackを`replica_XX/error.log`へ書き出します。
架橋構築で失敗した場合は、`crosslinking.py`のエラー行と標準エラー出力の末尾も含めます。
ランチャーは原因の1行を`jobs/status.json`の`error`へ記録し、今回の試行分のログ末尾300行を
`jobs/replica_XX.error.log`へ保存します。SIGKILLなどのシグナルで強制終了された場合は
シグナル名を記録します（SIGKILLはメモリ不足時のOOM killerでも発生します）。
構築結果は`replica_XX/crosslinking/results/epoxy/`、
その入力・到達反応率・構造ハッシュは`replica_XX/network_origin.json`に記録します。

同じprojectへの二重投入、および同じレプリカの二重実行はロックで拒否します。
失敗したレプリカがあっても他のレプリカは継続し、全ジョブ終了時にレポートとCSVを集計します。
`completed`は計算処理の正常終了を示し、物性値の収束判定とは別です。
再度`submit`すると完了済みレプリカを省略し、未完了レプリカを再開します。
完了済みの架橋構築はハッシュを検査して再利用します。
架橋構築そのものには途中再開機能がないため、中断した構築は同じseedで最初からやり直します。
その際は当該レプリカの未完了構築出力を作り直し、局所Sageテンプレートのキャッシュは再利用します。
最大cycleまでに反応率の目標へ届かなくても、その実現構造で後続評価を行います。
到達反応率と未達フラグをレポートに示し、未達レプリカはチェック通過値の集計に含めません。
別projectへ同じGPUを指定した場合のGPU共有は制限しないため、割り当てはproject間で調整してください。

単一レプリカを手動で実行することもできます。

```bash
CUDA_VISIBLE_DEVICES=0 python workflows/network_properties/network_properties.py run \
  workflows/network_properties/output/epoxy_gpu --replica 1
```

通常の`run`は従来どおり全レプリカを順次実行します。
集計ファイルにもロックをかけ、並列ジョブや手動`report`との書き込み競合を防ぎます。
`--smoke --replicates 2`のように、動作確認でも複数レプリカを明示できます。

## 熱履歴と300 K物性

| 段階 | 既定条件・出力 |
|---|---|
| 初期緩和と加熱 | 最小化、300 Kから100 K刻みで加熱、各100 ps |
| 高温緩和 | 800 K、NPT、1 ns |
| 冷却 | **800→100 K**、20 K刻み、各500 ps、平均40 K/ns |
| 冷却解析 | 各温度plateauの後半を使用。密度の双曲線フィットからTg、線・体積CTE |
| 温度復帰 | 100→300 K、20 K刻み、各100 ps |
| 300 K再緩和 | NPT、1 ns |
| 乾燥NPT観測 | 300 K、1 atm、5 ns、1 ps間隔 |
| 乾燥NVT観測 | 冷却・復帰後の固定セルで1 ns再緩和、5 ns観測 |
| 独立seed | 既定では架橋構築から独立した4レプリカ。保存構造入力時は共通トポロジー |

冷却は温度plateauを使う段階冷却です。800 Kにも観測用plateauを置くため、
冷却区間の合計は既定18 nsです（800→100 Kの温度差に対応する分は17.5 ns）。
時間は刻み幅とタイムステップに丸め、各段階の`done.json`に実時間を記録します。
既定タイムステップは0.5 fs。体積制御はwheelのconstraint group座標に対する
周期架橋対応Monte Carlo更新です。一般的な分子単位barostatへの置き換えはしません。

**100 Kまで冷却した後に300 Kへ戻した構造**を、弾性・溶媒評価に使います。
冷却途中の300 K構造も別名で保存します。再加熱でTgを再評価する設定ではありません。
800 Kで化学的に分解する場合でも、固定トポロジーMDには分解機構は含まれません。

RadonPyと同じ熱力学的揺らぎ関係を使う出力は次の通りです。

| 物性 | 関係・単位 |
|---|---|
| 密度 | 質量 / 体積、g/cm³ |
| 定圧比熱Cp | Var(H) / (kB T² m)、J/(kg K) |
| 等温圧縮率κT | Var(V) / (kB T ⟨V⟩)、1/Pa |
| 等温bulk modulus KT | 1/κT、GPa |
| 体積CTE αV | Cov(V,H) / (kB T² ⟨V⟩)、1/K |
| 線CTE αL | αV / 3（等方性を仮定）、1/K |
| 定積比熱Cv | Cp − T αV² / (ρ κT)、J/(kg K) |
| 断熱圧縮率κS | κT Cv / Cp、1/Pa |
| 断熱bulk modulus KS | 1/κS、GPa |
| bulk音速 | √(KS/ρ)、m/s。固体の縦波音速とは異なる |
| 直接Cv | NVTでのVar(E) / (kB T² m)、J/(kg K) |
| 固定電荷モデルの誘電応答 | NVTの連続的電荷変位から誘電テンソル、等方平均 |
| 幾何補助量 | wheelのpoint-probe空隙率など。溶媒アクセス可能体積とは異なる |
| ネットワーク形状 | 連結成分、最大成分の質量割合、cycle rank、周期セルを貫く結合cycleの有無 |

Hは`U + K + Pext V`です。比熱は古典・拘束付きSageモデルの値で、量子補正は行いません。
5ブロックの誤差、体積・エネルギーの有効サンプル数と前後半drift、
replica間の標本SDとSEMを出力します。1 replicaだけの場合、replica間誤差はnullです。
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

## 弾性

弾性は300 K構造から6成分×正負0.5%のひずみへ分岐します。
各分岐で100 ps NVT緩和、500 ps観測を行い、全ポテンシャルの微小セル変形差分
（PMEも含む）と運動項から応力を計算します。Voigt順序はxx, yy, zz, yz, xz, xy、
shearは工学ひずみです。

出力は6×6応力–ひずみ接線行列、ブロック誤差、対称性のずれ、
Voigt/Reuss/Hillのbulk/shear modulus、HillのYoung率・Poisson比です。
行列に非正の固有値があれば等方弾性率を確定しません。
仮想ひずみに拘束のvirialを欠落させないため、弾性計算では結合拘束を外し、
自然質量・0.25 fsを使います。NPTの拘束付きKTとの違いにはこのモデル差も含みます。

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
crosslinking/                       # 既定モード：コピーした反応・開始点ルール
input.structure.json               # --networkモードのみ：共通入力構造
summary.json / report.md
replica_properties.csv              # 各物性を行、各replicaの値・確認状態を列にした一覧
property_summary.csv                # 平均・標本SD・SEM・有効件数・チェック通過件数
jobs/                               # submit時：投入計画・状態・レプリカごとの実行ログ
replica_01/
  execution.json                    # レプリカ番号・seed・実行時に見えるGPU
  network_origin.json               # 既定モード：独立構築の反応率・原子数・構造ハッシュ
  crosslinking/results/epoxy/        # 既定モード：このレプリカ専用の架橋構造・履歴・設定
  cooling.csv
  cooling_300K.structure.json
  relaxed_100K.structure.json
  relaxed_300K.structure.json
  dry_result.json
  00_... / ... / 06_production_300K/  # samples, state, structure, done marker
  07_nvt_300K/                      # direct Cv and fixed-charge dielectric response
  mechanical/result.json
  solvents/water/
    packed/
    infinite_dilution_free_energy/  # windows, overlap, checkpoints, result
    diffusion/                     # MSD, diagnostics, OVITO frames, final structure
    result.json
```

同じ`run`コマンドで再開します。熱履歴は完了済みの温度・段階のstateから戻し、
中断した段階を初めから実行します。弾性は完了した正負分岐を
再利用します。溶媒MD・自由エネルギーはwheelのcheckpointから再開します。
熱履歴の再開はOpenMMの乱数状態まで一致する厳密再現ではありません。
`--network`で共通構造を入力した場合は、異なる架橋トポロジーによるばらつきを含みません。

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
OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python workflows/network_properties/network_properties.py new \
  workflows/network_properties/output/smoke \
  --network workflows/crosslinking/results/smoke/final.structure.json \
  --solvent water=O --solvent-molecules 1 --platform CPU --smoke
OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python workflows/network_properties/network_properties.py run \
  workflows/network_properties/output/smoke

python -m unittest discover -s tests -p test_network_properties.py -v
PAVES_NETWORK_INTEGRATION=1 OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m unittest discover -s tests -p test_network_properties.py -v
```

**`--smoke`は全段階を極端に短縮する接続テストです。値に物理的意味はありません。**
reportにも明記し、収束済みにはしません。

## 参照

- [RadonPy公式：対象物性と評価手順](https://github.com/RadonPy/RadonPy)
- [Hayashi et al., RadonPy, npj Comput. Mater. 8, 222 (2022)](https://doi.org/10.1038/s41524-022-00906-4)

RadonPyを依存として実行するのではなく、物性範囲と揺らぎ式を参照し、
PAVESが保存したSageネットワーク上で評価します。
