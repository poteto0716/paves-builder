# 保存電荷を保持した溶媒packと吸水MD

乾燥架橋の800→100 K熱緩和、300 K物性、溶媒和自由エネルギー・拡散係数を
一連で評価する場合は[架橋ネットワーク物性ワークフロー](../network_properties/README.md)を使用します。

[架橋構築](../crosslinking/README.md)の`final.structure.json`を読み込み、
樹脂の保存電荷・結合・原子IDを保持して溶媒を追加します。水が既定です。
溶媒1分子の座標、電荷、結合・角度・二面角・LJパラメータはPAVESのSage 2.3.0で
自動割り当てし、指定数を周期セルへpackします。

環境とwheelは架橋構築と共通です。リポジトリのルートで実行します。

```bash
python workflows/sorption/sorption.py pack \
  --network workflows/crosslinking/results/epoxy/final.structure.json \
  --smiles O --name water --molecules 100 --cell-scale 1.1 \
  --output workflows/sorption/results/water_100/packed

python workflows/sorption/sorption.py simulate \
  --packed workflows/sorption/results/water_100/packed \
  --output workflows/sorption/results/water_100/md --platform CUDA
```

`--smiles CO`ならメタノール、`--smiles 'CC(=O)C'`ならアセトン等にも変更できます。
Sageが扱える1つの連結した中性分子が対象です。イオンや混合SMILESは対象外です。
水もSageのモデルで、TIP3P等を選ぶ設定はまだありません。

packの既定`--cell-scale 1`では樹脂座標とセルを保持します。空間を広げるときは
両者を同じ倍率で拡張します。溶媒分子は周期境界をまたげます。
反対側の周期像も含めて重なりを検査し、指定数を配置できなければ停止します。
既存のpack出力は上書きしません。

MDは分子数を固定し、最小化と短い準備緩和の後、300 K・1 atmで100 ps NVT、
1 ns NPT緩和、5 ns NPT観測を行います。周期架橋に対応した体積更新を使います。
時間を変更する場合は`--nvt-ps`、`--npt-equilibration-ps`、`--production-ps`を指定します。
短いCPU確認は次のように実行できます。

```bash
OPENMM_CPU_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python workflows/sorption/sorption.py simulate \
  --packed workflows/sorption/results/water_100/packed \
  --output workflows/sorption/results/water_100/smoke --platform CPU \
  --nvt-ps 0.01 --npt-equilibration-ps 0.01 --production-ps 0.02 \
  --sample-interval-ps 0.01 --frame-interval-ps 0.01 --minimization-iterations 100
```

途中停止したMDは、同じコマンド・同じ条件に`--resume`を加えてcheckpointから再開します。
完了済み計算の再開は保存した結果を返します。条件が異なる再開は拒否します。

pack出力は`packed.structure.json`、`packed.system.xml`、`packed.settings.json`と
樹脂・溶媒の再利用テンプレートです。MD出力には時系列、体積・密度・充填率、MSD、
OVITO frames、最終構造、checkpoint、収束診断を含みます。構造JSONの座標とセルはÅ、
APIの`positions_nm`と`box_nm`はnmです。

**固定した水の量は平衡吸水率ではありません。** 質量基準の充填率は
`100 × 水の質量 / 乾燥樹脂の質量`です。今回は指定濃度での構造・膨潤・動力学を
計算するワークフローです。平衡吸水率には液体基準との化学ポテンシャル比較、
各自由エネルギー窓の収束、独立初期条件等の追加検証が必要です。
wheelにはアルケミカル自由エネルギーと比較解析のAPIも収録していますが、
このCLIは平衡吸水率を自動で確定しません。

Python APIだけでも同じ準備ができます。

```python
from polypaves import Solvent, pack_network
packed = pack_network("final.structure.json", Solvent("O", name="water"), 100, cell_scale=1.1)
packed.save("water_100")
```
