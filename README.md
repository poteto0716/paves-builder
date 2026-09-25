# PAVES

**Polymer Automation for Virtual Evaluation and Simulation**

*Pave the way from polymer chemistry to simulation.*

PAVES は、モノマーの SMILES と組成指定から分子動力学用の分子系を構築する
Python パッケージです。Python から C++ エンジンを直接呼び出し、座標・型・電荷・
結合・組成などを参照できる `System` を返します。正式名称は **PAVES**、Python パッケージ名と
コマンド名は **polypaves** です（`paves` は PyPI で別のパッケージが使っています）。
0.2.x までの名前は polyse でした。

この公開リポジトリは Linux 用バイナリ wheel、Python API の使用例、ドキュメント、
力場データを配布します。エンジンの C/C++ ソース、ヘッダー、オブジェクト、ビルド
ツリーは含みません。利用条件は [LICENSE](LICENSE) を確認してください。

## 対応環境とインストール

収録 wheel の検証環境は **Ubuntu 24.04 x86_64 / CPython 3.13** です。
Python ABI とプラットフォームが一致しない環境にはインストールできません。

```bash
git clone <repository-url> polypaves
cd polypaves
python3 -m pip install dist/polypaves-0.3.0-cp313-cp313-linux_x86_64.whl
python3 -c 'import polypaves; print(polypaves.__version__)'
sha256sum -c SHA256SUMS
```

実行時には Ubuntu 24.04 の `libc6`、`libstdc++6`、`libgcc-s1` が必要です。
コンパイラ、CMake、開発用ソースは不要です。

## 密着ワークフロー（高分子 / シリカ）

任意の高分子について、バルク → 表面 → シリカ上への加圧密着 → 緩和 → 界面エネルギー →
引張までを OpenMM で一続きに実行し、LAMMPS 形式でも出力するワークフローを
[`workflows/adhesion/`](workflows/adhesion/README.md) に収録しています。

```bash
conda env create -f workflows/adhesion/environment.yml && conda activate polypaves-adhesion
pip install dist/polypaves-0.3.0-cp313-cp313-linux_x86_64.whl
cd workflows/adhesion
python adhesion.py new pmma --monomer '*CC(C)(C(=O)OC)*' --dp 200 --chains 20
python adhesion.py run projects/pmma            # 中断しても同じコマンドで続きから
```

## 最初の分子系

```python
import polypaves

system = polypaves.build(
    monomer="*CC(C)(C(=O)OC)*",  # PMMA
    forcefield="pcff",
    chains=20,
    dp=100,
    density=1.18,
    temperature=413,
)

print(system.n_atoms)
system.write_lammps("pmma_system")
```

`dp` は1鎖の繰り返し単位数、`chains` は鎖数です。`write_lammps()` は指定
ディレクトリへ `system.data`、`system.in.styles`、`system.identity` を出力します。

## オブジェクトを組み合わせて pack する

`Polymer`、`Copolymer`、`Solvent`、`Slab` を組み合わせられます。

```python
import polypaves

pmma = polypaves.Polymer("*CC(C)(C(=O)OC)*", dp=20, name="pmma")
toluene = polypaves.Solvent("Cc1ccccc1", name="toluene")

# 成分ごとの個数を指定
system = polypaves.pack(
    [pmma, toluene],
    counts={"pmma": 4, "toluene": 100},
    forcefield="pcff",
    density=1.0,
)

# 全原子数の目標と weight 比を指定
system = polypaves.pack(
    [pmma, toluene],
    total_atoms=10_000,
    weight_fractions={"pmma": 0.7, "toluene": 0.3},
    forcefield="pcff",
    density=1.0,
)

print(system.composition)  # 実現した個数、mol比、weight比
```

`System` から `positions`、`atom_types`、`masses`、`charges`、`bonds`、
`angles`、`dihedrals`、`impropers`、`identity`、`box`、`composition`、
`report` を参照できます。詳しくは [Python API](docs/python_api.md) を参照してください。

## examples

[examples](examples/README.md) は `projects/paves/examples` の検証例を Python API
へ移植したものです。収録した実行入力はすべて `.py` で、独自入力ファイルや
CLI サブプロセスを使いません。

```bash
python examples/01_pcff_polystyrene/pmma_dp10/build.py
```

出力は各例の `output/` に生成され、Git 管理対象外です。力場の出典と収録データは
[外部データ](docs/external_data.md)、検証範囲は
[バイナリ検証](docs/validation.md) を参照してください。

## バイナリ配布と秘匿性

wheel 内の計算エンジンと高水準 API 実装は strip 済みネイティブ拡張です。公開側には
元の実装ソースを収録していません。一方で、利用者のコンピューターで動くバイナリの
解析を技術だけで完全に不可能にすることはできません。シンボル・デバッグ情報・私有
パスを除去し、[LICENSE](LICENSE) で逆コンパイル、逆アセンブル、実装復元を禁止する
構成です。詳細は [配布方針](docs/binary_distribution.md) を参照してください。
