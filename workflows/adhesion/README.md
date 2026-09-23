# 高分子 / シリカ 密着ワークフロー（PC-IFF, OpenMM）

任意の高分子について、バルク作成 → 表面作成 → シリカ上への設置 → 加圧密着 → 冷却・緩和
→ **界面エネルギー評価** → **引張（剥離）** までを一続きで実行します。構造と力場は
polyse（`dist/` の wheel）、分子動力学は OpenMM で計算し、各段階の構造は LAMMPS でも
そのまま実行できる形式で出力します。

## 準備（初回のみ）

```bash
conda env create -f workflows/adhesion/environment.yml
conda activate polyse-adhesion
pip install dist/polyse-0.2.0-cp313-cp313-linux_x86_64.whl
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
| 09 | pull | 上位 25%（正規化 z ≥ 0.75）原子群の重心をバネ 100 kcal/mol/Å² で 10 m/s（0.1 Å/ps）引張、150万 step。力を毎 step 出力、trajectory 保存 | 0.25 fs / なし |

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
| `05_compress/trajectory.dcd`, `09_pull/trajectory.dcd` + `topology.pdb` | trajectory（VMD, OVITO, MDTraj） |
| `09_pull/pull_force.csv` | 毎 step の基準点・重心・伸び・力（kcal/mol/Å と nN、+ は上向き） |
| `07_relax/fixed_atoms.lammps_ids.txt` | 固定原子の LAMMPS id（LAMMPS では `fix setforce 0 0 0` で固定） |

OpenMM の力場は `runs/01_build/polymer.openmm_system.xml`（高分子のみ）と
`runs/04_assemble/out/interface.openmm_system.xml`（界面系）です。

## LAMMPS ⇔ OpenMM 変換

```bash
python openmm_lammps.py to-lammps --state runs/07_relax/state.xml \
    --template runs/04_assemble/out/interface.data --styles runs/04_assemble/out/interface.in.styles --out relaxed.data
python openmm_lammps.py to-openmm --data relaxed.data \
    --system runs/04_assemble/out/interface.openmm_system.xml --out state.xml
```

## 注意

* OpenMM に slab 補正はないため、「z 方向の周期を外す」は十分な真空を入れた周期セルで
  扱っています。厳密な slab Ewald が必要なら、LAMMPS 出力を `boundary p p f` と
  `kspace_modify slab 3.0` で計算し直してください。
* 力場ファイル（PCFF, IFF）は `external/` の第三者データです。利用条件は `external/` の
  README を確認してください。
