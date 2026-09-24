# 高分子 / シリカ 密着ワークフロー（PC-IFF, OpenMM）

任意の高分子について、バルク作成 → 表面作成 → シリカ上への設置 → 加圧密着 → 冷却・緩和
→ **界面エネルギー評価** → **引張（剥離）** までを一続きで実行します。構造と力場は
polyse（`dist/` の wheel）、分子動力学は OpenMM で計算し、各段階の構造は LAMMPS でも
そのまま実行できる形式で出力します。

## 準備（初回のみ）

```bash
conda env create -f workflows/adhesion/environment.yml
conda activate polyse-adhesion
pip install dist/polyse-0.2.1-cp313-cp313-linux_x86_64.whl
```

GPU（CUDA または OpenCL）があれば自動で使います。無ければ CPU で動きますが、
6 万原子規模では非常に遅くなります。

## 新しい高分子系を投入する

```bash
cd workflows/adhesion

# 1. プロジェクトを作る（例: PMMA, DP 200 × 20 本 = 約 6 万原子, シリカ 3 × 2）
python adhesion.py new pmma --monomer '*CC(C)(C(=O)OC)*' --dp 200 --chains 20

# 2. 実行（途中で止まっても、もう一度同じコマンドで続きから再開）
python adhesion.py run projects/pmma

# 3. 進み具合と結果
python adhesion.py status projects/pmma
python adhesion.py report projects/pmma      # projects/pmma/results.html（3D表示とグラフ）
```

まず小さな系で通しを確認するのが安全です。`--test` を付けると全段階のステップ数が 1/100 になります。

```bash
python adhesion.py new ps_test --monomer '*CC(c1ccccc1)*' --dp 40 --chains 8 --supercell 1 1 --test
python adhesion.py run projects/ps_test
```

| `new` のオプション | 意味 |
|---|---|
| `--monomer SMILES` | 繰り返し単位（両端 `*`） |
| `--terminator SMILES` | 末端基（既定 `*C` = メチル） |
| `--dp N`, `--chains N` | 重合度と鎖数 |
| `--polyse-lines FILE` | 共重合体などを polyse の化学指定行で書く場合（`--monomer/--dp/--chains` の代わり） |
| `--supercell NX NY` | シリカのスーパーセル（既定 3 2 = 120.9 × 82.9 Å） |
| `--density` | ビルド密度 g/cm³（既定 1.0。バルク NPT で平衡密度へ） |
| `--forcefield` | `pcff-iff-mod`（既定。下の「修正 PCFF-IFF」参照）または `pcff-iff`（公開版そのまま） |
| `--test` | 全段階 1/100 ステップ |

すべての計算条件は `projects/<名前>/project.json` に書かれます。温度、ステップ数、
圧力、バネ定数、引張速度などを変える場合はこのファイルを編集してから `run` します。

## 途中から開始する

```bash
python adhesion.py run projects/pmma                    # 未完了の最初の段階から最後まで
python adhesion.py run projects/pmma --from compress    # compress 以降をやり直す
python adhesion.py run projects/pmma --to relax         # relax まででいったん止める
python adhesion.py run projects/pmma --only interface   # 1 段階だけ
```

各段階は前の段階の `runs/<NN_段階>/state.xml` から始まります。完了した段階には
`summary.json` があり、`status` で確認できます。前の段階が完了していなければ、
理由を示して止まります。

## 段階とプロトコル（既定値）

| # | 段階 | 内容 | dt / SHAKE |
|---|---|---|---|
| 01 | build | polyse で高分子メルトを作成。セル x, y = シリカスーパーセル | — |
| 02 | bulk | 変位制限付き予備緩和 → z 方向のみ可動の NPT (1 atm): 550 K 20万, 550→300 K 20万, 300 K 10万 step | 1 fs / あり |
| 03 | surface | 分子を切らずに z 方向へ真空 60 Å、NVT 550→300 K 20万 step | 1 fs / あり |
| 04 | assemble | シリカ上 3.0 Å に設置、上に真空 60 Å。polyse が結合系を作成 | — |
| 05 | compress | シリカ最下層 3.0 Å 固定、壁で 200 MPa、550 K 60万 step（trajectory 保存） | 1 fs / あり |
| 06 | cool | 壁を解除、550→300 K 20万 step | 0.25 fs / なし |
| 07 | relax | 300 K 50万 step | 0.25 fs / なし |
| 08 | interface | 300 K 1万 step、1000 step ごとに界面エネルギーを評価して平均 | 0.25 fs / なし |
| 09 | pull | 上位 25%（正規化 z ≥ 0.75）原子群の重心をバネ 100 kcal/mol/Å² で 10 m/s（0.1 Å/ps）引張、150万 step。力を 100 step ごとに出力（`pull.force_every`）、trajectory 保存。高分子の z 方向には熱浴をかけない（下記） | 0.25 fs / なし |

### 引張（09）の熱浴

Langevin 熱浴は実験室系で速度を 0 に引き戻すので、全原子にかけると、引き上げられて動く
高分子全体に −Mγv の抗力がかかり、バネはそれも支えることになります（PMMA DP 200 × 20 本、
γ = 1/ps、10 m/s で約 7 nN。剥がれたあとも力が 0 に戻らない原因）。そこで引張段階では
**高分子原子の z 成分だけ摩擦と乱数を外し**、x, y とシリカにだけ熱浴をかけます
（LAMMPS の `compute temp/partial 1 1 0` と同じ考え方）。z に入った熱は原子間力を通じて
x, y へ流れて抜けます。`project.json` の `pull.default_thermostat` を `"all"` にすると
従来どおり全方向にかけます。

検証（PMMA、剥離後の状態から 50 ps 継続）: 全方向 105.2 kcal/mol/Å（7.3 nN） →
z 除外 29.1 ± 2.7 kcal/mol/Å（2.0 nN）。残りはまだ界面に残る鎖を引き抜く実際の力です
（37.5 Å 引いた時点で 20 本中 7 本がシリカに接触）。完全に剥がれるまで見るなら
`pull.steps` を増やしてください。

### 界面エネルギー（08）

各スナップショットで、全系・高分子のみ・シリカのみの非結合エネルギー（LJ + Coulomb）を
**同じ座標・同じセル**で倍精度計算し、

E_int = E(全系) − E(高分子) − E(シリカ)、 γ = E_int / (Lx·Ly)

を求めて平均します（標準誤差つき）。結合項は界面をまたがないので厳密に打ち消し合います。
結果は `runs/08_interface/interface_energy.csv`（各スナップショット）と `summary.json`
（平均 ± 標準誤差、vdW / Coulomb 内訳、付着仕事 = −γ）にあります。

参考（PMMA DP 200 × 20 本、シリカ 3 × 2、10 スナップショット）:
γ = −127.2 ± 0.3 mJ/m²（E_int = −1,835 ± 5 kcal/mol; vdW −1,080, Coulomb −755）。

## 出力

各段階のディレクトリ `projects/<名前>/runs/<NN_段階>/`:

| ファイル | 内容 |
|---|---|
| `state.xml` | OpenMM State（座標・速度・セル） |
| `final.data`, `in.styles` | LAMMPS（`include in.styles` の後に fix と run を追加）。画像フラグ付き、速度込み |
| `final.pdb` | 表示用 |
| `log.csv`, `summary.json` | 経過と結果 |
| `05_compress/trajectory.dcd`, `09_pull/trajectory.dcd` + `topology.pdb` | trajectory（VMD, MDTraj。OVITO は下記） |
| `09_pull/pull_force.csv` | `pull.force_every` step ごと（既定 100。基準点は毎 step 動かす）の基準点・重心・伸び・力（kcal/mol/Å と nN、+ は上向き） |
| `07_relax/fixed_atoms.lammps_ids.txt` | 固定原子の LAMMPS id（LAMMPS では `fix setforce 0 0 0` で固定） |

OpenMM の力場は `runs/01_build/polymer.openmm_system.xml`（高分子のみ）と
`runs/04_assemble/out/interface.openmm_system.xml`（界面系）です。

### OVITO で見る

```bash
python dcd_to_ovito.py projects/<名前>/runs/09_pull      # compress なら runs/05_compress
```

`trajectory_ovito.dump`（`id x y z` のみ、セルは `final.data` と同じ 0〜L）ができます。
OVITO で `final.data` を LAMMPS data（atom style `full`）として開き、**Load trajectory**
モディファイアで `trajectory_ovito.dump` を読み込むと、`final.data` の原子タイプ・電荷・
結合を保ったまま座標だけがフレームごとに変わります。MDTraj の `save_lammpstrj` の出力は
全原子を type 1 と書くため、Load trajectory でタイプが上書きされて色分けできなくなります。

## LAMMPS ⇔ OpenMM 変換

```bash
python openmm_lammps.py to-lammps --state runs/07_relax/state.xml \
    --template runs/04_assemble/out/interface.data --styles runs/04_assemble/out/interface.in.styles --out relaxed.data
python openmm_lammps.py to-openmm --data relaxed.data \
    --system runs/04_assemble/out/interface.openmm_system.xml --out state.xml
```

## 修正 PCFF-IFF（既定、`--forcefield pcff-iff-mod`）

既定の力場は `examples/forcefields/pcff_iff_long_bulk_mod.ff` です。公開版の PCFF-IFF を
そのまま使うには `new --forcefield pcff-iff` とします。

PCFF の Si–O–Si（`sio-osi-sio`）は θ0 = 157°、180° でも 0.42 kcal/mol しか上がらない
非常に柔らかい角度で、熱運動で 179.9° 以上まで開きます。この角をまたぐ二面角は直線で
定義できず、力が 1/sin θ で発散します（179.9° 以上で O に最大 7,400 kcal/mol/Å）。
1 fs + SHAKE では、直線を横切るたびの大きな力でエネルギーが保存されず、ときどき拘束が
破綻して止まります。0.25 fs で止まらないのは、特異点付近を細かく追えるからです。

`pcff_iff_long_bulk_mod.ff` は、θ0 ≥ 150° の角をもつ二面角だけについて、φ に依存する項
（ねじれ、middle/end-bond-torsion、angle-torsion、angle-angle-torsion）に
S(θ1)·S(θ2) を掛けます。S は 175° までは厳密に 1、175° → 180° で
1 − x³(10 − 15x + 6x²) により滑らかに 0 になります。係数は変えず、角度項（シロキサンの
柔軟性）もそのままです。

メチルフェニルシロキサン（8,288 原子、Si–O–Si 288 個）、550 K、1 fs + SHAKE での比較:

| | PCFF-IFF | 修正 |
|---|---|---|
| 179.9° 以上での O の最大力 | 7,397 kcal/mol/Å | 112 kcal/mol/Å |
| NVE 20 ps の全エネルギー（ドリフト / RMS） | +142 kcal/mol/ps / 844 kcal/mol | −0.29 kcal/mol/ps / 1.9 kcal/mol |
| Si–O–Si の平均 ± 標準偏差 | 150.90 ± 9.66° | 150.83 ± 9.75° |
| 175° を超える割合 | 0.87 % | 1.16 % |
| 800 K、1 fs + SHAKE、50 ps | 12 本中 1 本 NaN | 8 本中 0 本 |

175° 未満の構造では、エネルギーも力も元の PCFF-IFF と一致します（相対 2×10⁻¹⁶）。
θ0 ≥ 150° の角がない系（PMMA と IFF シリカなど）では、書き出される System は同一です。
LAMMPS の `dihedral_style class2` にはこの切り替えがないため、LAMMPS 出力は元の
PCFF-IFF のままです（`in.styles` に注記が入ります）。

## 注意

* OpenMM に slab 補正はないため、「z 方向の周期を外す」は十分な真空を入れた周期セルで
  扱っています。厳密な slab Ewald が必要なら、LAMMPS 出力を `boundary p p f` と
  `kspace_modify slab 3.0` で計算し直してください。
* 力場ファイル（PCFF, IFF）は `external/` の第三者データです。利用条件は `external/` の
  README を確認してください。
