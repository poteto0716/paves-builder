# Epoxy-only exampleのYAMLを書く・変更する

このexampleは、未反応epoxideの一部を中性OH開始点へ変換し、そのOHから
`OH* + epoxide → elongated OH*` を繰り返すものです。
`*` は `reaction_state: cationic_activity` というReaction Engine専用metadataです。
formal charge、radical、OpenFF partial chargeを意味しません。

**現在のYAMLですべてのエポキシの計算成功を保証するものではありません。**
検索SMARTSへの一致、意図した開環位置、初期構造と生成物のSage割り当ては、それぞれ確認が必要です。
対応範囲と変更時の確認方法を後半にまとめています。

## 1. どのファイルに何を書くか

| 場所 | 書く内容 | このexampleの設定 |
|---|---|---|
| [activation.yaml](activation.yaml) | MD前に開始点を作る方法、対象官能基、開始点数／割合、追加原子、結合編集、state付与 | epoxideを1本開環し、CとOをHでcapしてOH*を作る |
| [reactions.yaml](reactions.yaml) | MD中に繰り返す反応、反応前後のSMARTS、state条件、距離条件、反応上限 | OH*が次のepoxideを開環し、activityを新OHへ移す |
| [crosslinking.py](crosslinking.py) | モノマーSMILES・総原子数・密度、Sageバージョン、初期NPT・cycle NVT/NPT、CUDA、出力先 | 3官能モノマー、20000原子超、目標反応率99.5%、最大5000 cycle、Sage 2.3.0、CUDA/mixed |

library exampleは、実行時の**カレントディレクトリ**の `activation.yaml`／`reactions.yaml` を読みます。
`--project 名前` で `results/名前/` に保存し、両YAMLもコピーします。
既存projectの再実行には `--force` が必要です。指定projectの内容を削除して作り直します。

例えば別のエポキシへ変更する場合、まずPythonの `p.Solvent(...)` のSMILESを変更します。
反応させる官能基と反応前後の局所構造が同じなら、YAMLを変更せずに使える可能性があります。
サイズは `--total-atoms`、密度は `--density`、SMILESは `--smiles`／`--smiles-file` で指定します。
総原子数はHと開始点作成で追加するHを含み、分子単位でabove/nearest/floor/exactを選びます。
初期緩和は `--equilibration-steps`／`--equilibration-timestep`、成長後のMDは
`--nvt-steps`／`--nvt-timestep` と `--npt-steps`／`--npt-timestep` で指定します。
初期NPT緩和の既定は300 K・1 atm、1 fs × 1,000,000 step（1 ns）です。
OVITO軌跡とcycle別構造は既定では保存せず、`--save-trajectory`で有効にできます。
有効時は各cycle終了時に1回保存します。反応履歴・集計JSON・最終構造はOFFでも保存します。
これらをこの2つのYAMLへ追加すると、未対応の設定キーとしてエラーになります。

なお、現在の短いlibrary exampleにはPython側の明示指定があります。
YAMLを編集する前に、後述の「PythonとYAMLの優先関係」を確認してください。

## 2. YAMLそのものの書式

```yaml
# #以降はコメント。ここはプログラムに読まれない。
activations:                         # キー: 値。コロンの後に空白を入れる
  - name: neutral_ring_opening_seed  # - はリストの1項目。ここでは1つのルール
    target_smarts: "[C;H2:1]1[O;+0;X2:2][C:3]1"
    required_atom_states:
      2: null                       # map 2にはstateがないことを要求
    selection: fraction
    fraction: 0.05                  # 数値。5%は0.05と書く
```

- インデントはスペースでそろえます。このexampleでは1階層2スペースです。タブを使いません。
- 同じ階層のキーは同じ位置から書きます。`name` と `target_smarts` は同じルールの項目です。
- 複数のルールを書く場合は `- name: ...` を同じインデントで追加します。ルール名は同じファイル内で一意にします。
- SMARTSは文字列として `"..."` または `'...'` で囲みます。`#`、`[]`、`:` などをYAMLの記号と区別できます。
- 真偽値は `true` / `false`、整数は `1`、小数は `0.05` と書きます。数値や真偽値は引用符で囲みません。
- `null` は「値なし」です。`"null"` と書くと文字列のstate名になり、意味が変わります。
- `[1, 2]` はリスト、`{2: null}` は辞書です。複数行のリスト／辞書でも同じ内容を表せます。
- `fraction: 5` や `fraction: 5%` は使えません。割合の範囲は0〜1です。
- 同じキーを1つの辞書内に重複させません。例えば同じルール内に `fraction` を2回書かないでください。

数値のatom mapキーは `1: ...` の形で書いて構いません。JSONの場合は `"1": ...` と書き、loaderが整数へ変換します。
読み込みはYAMLでは `yaml.safe_load` を使います。Pythonの式、関数、コードを書く場所ではありません。

## 3. SMARTSの数字とstable atom IDの違い

初期化の検索パターンは次です。

```text
[C;H2:1]1[O;+0;X2:2][C:3]1
```

| 記号 | 意味 |
|---|---|
| `C` | 非芳香族の炭素を検索する |
| `H2` | 結合しているHが2個（explicit Hも数える）。末端CH2だけに一致する |
| `;` | 同じ原子への条件のAND |
| `O` | 非芳香族の酸素を検索する |
| `+0` | formal chargeが0であることを要求する |
| `X2` | Hを含む結合相手の数が2であることを要求する |
| `:1`、`:2`、`:3` | このルール内で原子を対応付けるatom map番号 |
| 角括弧の外側の `1 ... 1` | 環を閉じる記号。これによりC–O–C三員環を指定する |
| `[#1:5]` | 元素番号1、つまり水素。そのatom map番号は5 |
| `-` | 単結合 |

**atom map番号はOVITOのparticle ID、RDKit index、OpenMM particle indexではありません。**
毎回のマッチで、例えば「map 2 → この環の酸素のstable_atom_id」という対応が作られます。
同じruleを別の環に適用すれば、map 2が指す実原子も変わります。
`activation.yaml` と `reactions.yaml` のatom mapは独立です。
初期化の酸素がmap 2だからといって、成長ruleのmap 2も酸素になるわけではありません。

原子の実体は `stable_atom_id` で保存され、履歴にも記録されます。
同じactivityが別の酸素へ移った場合、移動元と移動先のstable IDは別ですが、activity数は保存します。
このexampleのstateは文字列であり、独立した `activity_id` をこのYAMLで指定する機能はありません。

## 4. activation.yaml：初期OH*を作る

実際のファイルには各項目のコメントを入れてあります。初期化は次の編集です。

```text
before: C1–O2–C3–C1 のepoxide環
after:  H4–C1–C3–O2–H5
                          *

remove: C1–O2
add:    C1–H4、O2–H5
state:  O2にcationic_activityを付与
```

環の外の置換基・結合は保持します。追加Hは選んだ環あたり2個です。
既存原子のstable IDは保存し、追加Hには新しいIDを割り当てます。
原子追加はMD前の初期化で行い、その後の成長では原子を追加しません。
このH cappingは開始構造を準備するproxy操作で、PAG光化学の再現ではありません。

| キー | 書く値と意味 |
|---|---|
| `name` | 必須。履歴で区別するルール名 |
| `target_smarts` | 必須。初期化対象の単一connected SMARTS。対象原子に重複しないatom mapを付ける |
| `active_atom_map` | 必須。stateを付ける原子。この例では元epoxideのOである2 |
| `reaction_state_to_assign` | この例では `cationic_activity`。付与するgeneric state名 |
| `required_atom_states` | 初期化対象のstate条件。`2: null` はstateのない酸素だけを対象にする |
| `selection` | `count` / `fraction` / `random` / `explicit`。通常はcountかfractionを使う |
| `count` | 非負整数。選ぶ官能基の個数。fractionと同時に指定しない |
| `fraction` | 0〜1。選ぶ官能基の割合。countと同時に指定しない |
| `seed` | 非負整数。ランダム選択に使う乱数seed。Python側で上書きできる |
| `atoms_to_add` | 追加原子のリスト。新しいatom map、元素、初期座標の作り方を書く |
| `bonds_to_remove` / `bonds_to_add` | 切る／作る結合。`[map1, map2, bond_order]` で書く |
| `require_explicit_hydrogens` | この例では `true`。Hが実際の原子として存在する構造を要求する |
| `expected_total_charge_delta` | 選択対象1個あたりのtotal formal charge変化。この例では0 |
| `expected_atom_state_deltas` | 選択対象1個あたりのstate数変化。初期化では `cationic_activity: 1` |
| `metadata` | 原則として履歴・説明用。`representation` 等を書くだけで力場や化学処理が切り替わるわけではない |

count / fractionはいずれも、適格な対象からseedに基づいてランダムに選びます。
対称なSMARTSのマッチを同じ官能基としてまとめたうえで数えます。
このexampleでは**分子数ではなくepoxide環の個数**が分母です。
選択数は `floor(fraction × 適格なepoxide環数)` です。

```yaml
# 正確に1個を開始点にする：selection/fractionの2行を次に置き換える
selection: count
count: 1

# 割合で選ぶ：selection/countの2行を次に置き換える
selection: fraction
fraction: 0.05
```

この2つは代替設定です。同じルールに両方のブロックを書かないでください。
3官能×8分子なら24環なので、5%は1環になり、実際の開始割合は1/24です。
`count: 0` / `fraction: 0.0` では開始点ができず、この成長ruleは反応しません。
`fraction: 1.0` では全環をOH化し、成長に必要な未開環epoxideが残りません。
`count: 1` は全環ではなく1環だけです。

特定の原子を指定する場合は、酸素の**stable ID**で選びます。

```yaml
selection: explicit
selection_atom_map: 2
explicit_atom_ids: ["対象酸素のstable_atom_id"]
```

この場合はcount/fractionを削除します。文字列は説明用なので実際のIDへ置き換えます。
OVITOの表示番号をそのまま書かず、構造／履歴のstable IDとの対応を確認してください。

### 追加Hと位置指定

```yaml
atoms_to_add:
  - atom_map: 4
    element: H
    position:
      relative_to: 1
      towards: 2
      distance: 1.09
```

これは「既存map 1の位置からmap 2の方向へ1.09 Å進んだ位置に、新H4を置く」という指定です。
`towards` は結合を作る指示ではありません。結合は別の `bonds_to_add` で作ります。
周期境界がある場合、方向はminimum imageで決めます。追加直後の近接や歪みは最小化で緩和します。
`element: H` の代わりに `atomic_number: 1` も使えますが、両方を同時に書きません。

genericな初期化には次の編集キーも使えます。今回のepoxideでは必要ありません。

```yaml
bond_orders_to_change:
  - {atom_map_1: 1, atom_map_2: 3, order: 2, old_order: 1}
formal_charges_to_change:
  - {atom_map: 2, charge: 0}
atom_properties_to_change:
  - {atom_map: 1, property: user_tag, value: initialized}
```

ここは独立した**書式例**です。このepoxide ruleへそのまま追加しないでください。
`atoms_to_remove` は将来用で、現時点では空リスト `[]` のみです。原子削除は未実装です。
identity、atom type、partial chargeをcustom propertyで書き換えることはできません。
同じ原子へのreaction_state付与を `reaction_state_to_assign` とproperty編集で二重指定しません。

## 5. reactions.yaml：OHから次のepoxideへ伸ばす

このファイルの `reactions:` は成長ルールのリストです。
top-levelの `selection` / `seed` は選択方針と乱数設定を保持します。
この設定をrunnerへ渡す方法は次節に示します。

```yaml
selection: nearest
seed: 31415
reactions:
  - name: epoxy_OH_relay
    reactants:
      - "[O;+0;X2:1]-[#1:5]"
      - "[C;H2:2]1[O;+0;X2:3][C:4]1"
    products:
      - "[O;+0:1]-[C:2]-[C:4]-[O;+0:3]-[#1:5]"
    reactive_atom_maps: [1, 2]
    required_atom_states: {1: cationic_activity, 3: null}
    product_atom_states: {1: null, 3: cationic_activity}
    require_explicit_hydrogens: true
    distance_cutoff: 5.0
    allow_intramolecular: true
    minimum_graph_separation: 4
    max_reactions_per_cycle: 1
    expected_total_charge_delta: 0
    expected_atom_state_deltas: {cationic_activity: 0}
    metadata:
      conversion_reactant: 1
```

この抜粋も読み込み可能です。実ファイルにはproxyの説明metadataも記録しています。

| map | 反応前 | 反応後 |
|---|---|---|
| 1 | activityを持つ旧OH酸素 | 新しいエーテル酸素。stateを解除 |
| 2 | epoxideの攻撃される末端CH2炭素 | 旧OH酸素1と結合する炭素 |
| 3 | stateのないepoxide酸素 | H5を受け取る新OH酸素。activityを付与 |
| 4 | epoxideのもう一方の炭素 | C2とO3への結合を保つ |
| 5 | 旧OHに結合する既存H | 新OHに結合する同じH |

SMARTSの差分から、O1–H5 / C2–O3を切り、O1–C2 / O3–H5を作る編集を自動生成します。
`products` は局所断片です。巨大polymer全体を書き直す必要はなく、マッチ範囲外の構造は保持します。
成長ruleでは反応前後のmapped atom集合が同じでなければなりません。追加・削除はできません。
原子を追加する初期化の `atoms_to_add` を、成長のReactionRuleへ書くことはできません。

state条件はSMARTS条件と同時に判定されます。
このruleでは `OH* + neutral epoxide` だけが候補です。
普通のOH、既にstateのあるepoxideは候補になりません。
`required_atom_states` のキーを省略すると、その原子のstateは制約しないという意味になります。
`null` を書いて「stateなし」を要求することとは異なります。
`product_atom_states` の `null` はstateの解除です。

| 調整キー | 意味・注意点 |
|---|---|
| top-level `selection` | `nearest` は距離優先、`random` は乱数、`rate_weighted` はルールのrate_weightによる重み付け。1batchに1方針 |
| `reactive_atom_maps` | 現在座標で近接を見るアンカー。この2反応物ではO1–C2の距離。編集対象原子の一覧ではない |
| `distance_cutoff` | 近接判定の上限、単位Å。PBC minimum imageで比較する。既定値は省略時3.5、このYAMLは5.0 |
| `allow_intramolecular` | 現在同じconnected componentに属するサイト同士を許可するか。trueなら閉ループも許す |
| `allow_intermolecular` | 異なるconnected component同士を許可するか。省略時true。このexampleでもtrue |
| `minimum_graph_separation` | 同一componentの場合、アンカー間の最短結合経路をこの本数以上に制限。4はÅではなく結合4本 |
| `max_reactions_per_cycle` | このルールの1cycleあたり上限。正整数、またはnullでルール上限なし |
| `expected_total_charge_delta` | 1反応あたりのtotal formal charge差。0でもpartial charge再割り当ては行う |
| `expected_atom_state_deltas` | 1反応あたりのstate数差。0でactivity数保存を検証。違反時はtransactionをrollback |
| `metadata.conversion_reactant` | conversion対象のreactantリストindex。**0始まり**。1は2番目のepoxideを指定 |
| `maximum_conversion` | 必要なら0〜1でルールのconversion上限を指定。runner全体のtarget_conversionとは別 |
| `priority` / `rate_weight` | 複数ルールの優先度／重み。省略時0／1.0。化学反応速度定数の自動推定ではない |

`allow_intramolecular: false` の場合、一度同じcomponentになった分子の残り官能基同士も反応禁止になります。
ここでいう分子内／分子間は、初期モノマー単位ではなく**現在の結合グラフのcomponent**を基準にします。
同一atomを共有する候補は同じbatchに採用しません。stateが1個のこのruleではその共有も制限になります。

反応cutoffとOpenMM nonbonded cutoffは別の設定です。
探索範囲は平衡結合長や反応障壁ではありません。
値を小さくすると候補が減り、反応には接近を待つMDが必要になる場合があります。
crosslinking.pyはCLIの開始距離（既定5 Å）でYAMLの距離を上書きし、0件が2 cycle続くと
次cycleから1.2倍、最大10 Åまで自動拡大します。`--fixed-reaction-cutoff`で拡大を無効化できます。

conversionは、この例ではrunnerを作った時点（初期化後）の適格な未開環epoxideに対する成長時の消費割合です。
初期OH化されたepoxideはその分母に入りません。
24環中1環を初期化し18反応した場合、成長conversionは18/23、初期化込みの消費率は19/24です。
`mechanism`、`representation`、`forcefield_assumption` は説明・履歴用ですが、
`conversion_reactant` は実際に集計対象を指定する設定です。

## 6. PythonとYAMLの優先関係

**現在のcrosslinking.pyでは、YAMLの全項目が自動で実行設定になるわけではありません。**
開始点のCLI指定を省略するとYAMLの個数／割合を使います。reaction選択方針もYAMLから読みます。
cycleの反応数上限はCLIで上書きできます。

| YAML設定 | 現在のPython指定 | 有効になる値 |
|---|---|---|
| activationの `selection/fraction/count` | CLIのcount/fraction、またはYAML | CLI指定を優先。省略時はYAML（現在10%） |
| activationの `seed` | `ActivationInitializer(..., seed=seed)` | Python変数seed=31415 |
| reactionsのtop-level `selection` | `ReactiveMDRunner(..., selector=rules.selection)` | YAMLのselection（現在nearest） |
| reactionsのtop-level `seed` | `ReactiveMDRunner(..., seed=seed)` | Python変数seed=31415 |
| 各ruleの `max_reactions_per_cycle` | CLIの`--max-reactions-per-cycle` | 指定時はルール／全体上限を上書き。省略時はYAML |
| 各ruleの `distance_cutoff` | CLIの`--reaction-cutoff`と自動拡大条件 | 開始距離はCLI（既定5 Å）。既定では0件が2 cycle続くと1.2倍、最大10 Å |

反応後の電荷・力場の割り当て方式はYAMLではなくPython/CLIで指定します。
crosslinking.pyの既定は`--reaction-assignment local_templates`です。
更新5結合＋参照8結合の小フラグメントの電荷を再利用し、電荷総和を補正・検証します。
全体再計算は`--reaction-assignment whole`です。
[公開ワークフローの実行と出力](README.md)を参照してください。

YAMLの `fraction: 0.20` はCLIで個数／割合を指定しなければ有効になります。
CLIから割合を指定するなら `--activation-fraction 0.20` を使います。

YAML側で開始割合と選択方針を管理したい場合は、次のように該当部分を置き換えます。
`structure`、`forcefield`、`bridge`、`validate_initialization`、温度等はexampleで作ったものを使います。

```python
activation_rule = load_activation_definitions(definitions / 'activation.yaml')[0]
# seed引数を省略して、activation_rule.seedを使う
initializer = ActivationInitializer(structure, forcefield=forcefield)
# count/fractionを省略して、YAMLのselectionと値を使う
initializer.apply(activation_rule, validator=validate_initialization)

reaction_definitions = load_reaction_definitions(definitions / 'reactions.yaml')
runner = ReactiveMDRunner(
    structure=structure,
    forcefield=forcefield,
    reactions=reaction_definitions,
    selector=reaction_definitions.selection,
    seed=reaction_definitions.seed,
    temperature=temperature,
    timestep=nvt_timestep,
    platform=platform,     # exampleの既定はCUDA
    precision=precision,
)
```

`runner.run` で全体上限の `max_reactions_per_cycle=1` を外すか `None` にすれば、各ruleの上限を使います。
必要なreporter・履歴・最小化設定は元のexampleに従います。
初期化seed、packing seed、MD thermostat seedは別の用途です。
同じseedでも、異なるGPU/CPUや数値精度をまたいだ完全な軌跡一致を保証するものではありません。

## 7. すべてのエポキシにこのまま使えるか

確認は次の3段階です。

1. **検索範囲**：その官能基がC–O–C三員環SMARTSに一致するか。
2. **指定するgraph変換**：攻撃するC、切るC–O、Hの移動、state移動が目的のtopologyを表すか。
3. **計算可能性**：初期化後と各生成物にstock Sageを完全に割り当てられ、有限エネルギーでMDできるか。

小分子のRDKit検索では、現在のactivationとgrowthのepoxide SMARTSに次の結果を確認しました。
**この表の検索一致は、Sage/OpenMM実証済みという意味ではありません。**

| 構造例 | 検索一致する環数／分子 | Sage/OpenMM検証の範囲 |
|---|---:|---|
| ethylene oxide `C1CO1` | 1 | この表では検索のみ |
| propylene oxide `CC1CO1` | 1 | この表では検索のみ |
| methyl glycidyl ether `COCC1CO1` | 1 | この表では検索のみ |
| DGEBA `CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1` | 2 | この表では検索のみ。生成物のSage割り当ては別途確認 |
| TMP triglycidyl ether（このexampleのモノマー） | 3 | 8分子・1開始点・18成長反応、全stepでstock Sage割り当て／System再構築／有限エネルギー／CPU MDを検証済み |
| cyclohexene oxide `C1CCC2OC2C1` | 0 | 末端CH2を持たないので対象外。使うなら攻撃Cの条件を変える |
| tetrahydrofuran `C1CCOC1` | 0 | C–O–C三員環ではないので対象外 |
| aziridine `C1CN1` | 0 | 環原子がOではないので対象外 |

モノ・二・三官能という数そのものはReaction Engineの分岐条件ではありません。
1分子に条件を満たすepoxideが複数あれば、それぞれを官能基として検索します。
ただしネットワークが形成されるかは官能基数、成長の進み方、近接条件、分子内反応の許可などに依存します。
エポキシを含んでいるという理由だけで必ず架橋・gel化するとは限りません。

### 攻撃位置と立体化学

成長反応の攻撃C2と初期化の切断C1は、どちらも `[C;H2]`（末端CH2）に限定しています。
塩基性・求核的なSN2型開環（アミン硬化、アルコキシドなど）で、置換基の少ない炭素がほぼ選択的に攻撃される位置選択を表します。
そのためglycidyl etherの結合は常に `O–CH2–CH(OH)–` の向き（head-to-tail）になり、CH側への攻撃は起こりません。
初期化では切断したCH2がCH3になり、OHは第二級炭素CHに付きます。

酸性・カチオン条件のようにCH側への攻撃も混ぜたい場合は、両方のSMARTSを `[C:2]` / `[C:1]` に戻します。
その場合、異なる向きのマッチは候補の距離と選択方針で選ばれ、エネルギー的な優先性は入りません。
末端CH2を持たない内部epoxide（cyclohexene oxide等）は、`[C;H2]` のままでは検索されません。

また、現versionはproduct側の明示的な `@` 等のstereochemistry指定に対応していません。
反応時の新しい立体配置や反転を自動で決定する機能もありません。
検索一致したことだけで、立体化学まで目的どおりになるとは判断できません。

OH側のSMARTSも、現在は「neutral X2 oxygenとHの結合」という広い条件です。
水などを追加する拡張例で、炭素に結合するOHだけを検索したい場合には、
`[C,c]-[O;+0;X2:1]-[#1:5]` のようにread-onlyな炭素contextを加えられます。
mapのないcontext原子は生成物へ書かず、編集対象にもなりません。
今回のepoxide-only例では初期化で作ったOHだけがstateを持つため、追加のアルコール等は必要ありません。

### Sageの対応範囲

`+0` を書いて中性にしたことは、完全parameterizationの保証にはなりません。
PAVESの現OpenFF adapterはH/C/N/O/F/P/S/Cl/Br/Iを対象とし、例えばSiや金属はこの経路では対象外です。
この元素集合に収まっていても、対象の結合環境すべてにSage parameterがあるかは別途検証します。
モノマーだけでなく、初期OH化後・各開環生成物・分岐や閉ループを含む生成物を確認する必要があります。

このexampleはSage .offxmlを変更せず、GAFF2等へのfallbackもしません。
未割り当てparameter、valence、state数・formal charge差、System構築、非有限エネルギーが問題になると
reaction transactionは失敗し、graph／stateの更新をrollbackします。
説明用metadataの `forcefield_assumption` を変えても、欠けたparameterが補完されることはありません。

## 8. 読み込みだけを先に確認する

GPU計算前に、YAML書式・キー名・SMARTS compilationを確認できます。
このコードはpackingやMDを実行しません。Sageの完全割り当て確認は別途必要です。

```python
from pathlib import Path
from polypaves.reaction import load_activation_definitions, load_reaction_definitions

directory = Path.cwd()
activation = load_activation_definitions(directory / 'activation.yaml')[0]
reactions = load_reaction_definitions(directory / 'reactions.yaml')

print(activation.name, activation.selection, activation.count, activation.fraction)
print(reactions.selection, reactions.seed)
print(reactions[0].reaction_smarts)
print(reactions[0].graph_edit)  # 自動compileされた結合／state編集
```

新しいモノマーでは少数分子・少数反応から開始し、全生成物のtypingとenergyを確認してから系を大きくします。
MDのplatform指定はexampleでCUDAを使います。CPU確認をしたい場合だけ `--platform CPU` を明示します。
ここに記載した既存の18反応検証はCPUで実施し、CUDAでの実行は未検証です。

| 状況 | 最初に確認するもの |
|---|---|
| YAMLを読み込めない | インデント、SMARTSの引用符、数値型、未知キー、count/fractionの二重指定 |
| 初期化0件 | epoxide数に対するfloor計算、target SMARTS、required state、Python側の上書き |
| 成長候補0件 | activityを持つOHの存在、未反応epoxide、explicit H、距離cutoff、component条件、graph separation |
| 想定と違うcarbonに結合する | attack map 2のSMARTS条件。現在は `[C;H2]` で末端CH2だけを攻撃する |
| Sageで失敗する | 初期化後／生成物の具体的な結合環境と未割り当てterm。SMARTS一致とFF対応を区別する |
| activityが増えたように見える | OVITOの `ReactionActive` と履歴のstate数を確認。`ReactiveSite` は未反応epoxideも含む |

このneutral proxyは、指定したルールに基づくpolyether/network connectivityを作るための表現です。
実oxonium、ion pairing、PAG光化学、プロトンの実運動、反応障壁や重合速度は明示的に再現しません。
物理的な機構・速度の検証と、中性graph変換／FF更新が正常に動くことの検証は別です。
