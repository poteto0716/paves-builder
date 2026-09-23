# Python API

## `build()`

単一ポリマー用の簡便 API です。

```python
system = polyse.build(
    monomer="*CC*",
    dp=20,
    chains=4,
    terminator="*[H]",
    forcefield="examples/forcefields/opls_alkane.ff",
    density=0.8,
    temperature=413,
)
```

`chains` の代わりに `total_atoms` と `atom_count_mode="floor" | "nearest" |
"exact"` を指定できます。両端を変える場合は単一 SMILES、非対称なら
`terminator=(head, tail)`、キャップなしなら `None` を指定します。

## component objects と `pack()`

```python
polyse.Polymer(monomer, dp, name=None, terminator="*[H]")
polyse.Copolymer(monomers, dp=None, architecture="random_exact",
                 fractions=None, sequence=None, blocks=None,
                 name=None, terminator="*[H]")
polyse.Solvent(smiles, name=None)
polyse.Slab(path, name=None, format=None, fixed=True, preserve=False,
            exclude_radius=0.0, replicate=(1, 1, 1))
```

`pack()` のサイズ・組成指定は次のいずれかです。

| 指定 | 意味 |
|---|---|
| `counts={name: N}` | ポリマーは鎖数、溶媒は分子数 |
| `total_molecules=N` + `mole_fractions` | 合計分子数と鎖 mol 比 |
| `total_molecules=N` + `weight_fractions` | 合計分子数と質量比 |
| `total_atoms=N` + 比率 | 許容幅内で整数分子の組成を解く |

混合系の原子数許容幅は `tolerance=0.10` が既定です。分子や鎖を途中で切らない
ため、実現原子数は目標と一致しない場合があります。`system.requested` と
`system.report["packing"]` で要求値と丸め結果を確認してください。

### Copolymer

```python
copolymer = polyse.Copolymer(
    monomers={"M": "*CC(C)(C(=O)OC)*", "S": "*C(c1ccccc1)C*"},
    dp=40,
    fractions={"M": 0.7, "S": 0.3},
    architecture="random_exact",
    name="copolymer",
)
```

`random_exact`、`random_probabilistic`、`alternating` のほか、
`sequence=["M", "S", ...]` または `blocks=[("M", 10), ("S", 20)]` を
指定できます。

### Slab と充填領域

```python
surface = polyse.Slab("surface.mol2", name="surface", fixed=True, exclude_radius=2.0)
region = polyse.BoxRegion(origin=(0, 0, 25), lengths=(40, 40, 55))
system = polyse.pack(
    [surface, pmma, toluene],
    counts={"pmma": 4, "toluene": 100},
    forcefield="surface_and_polymer.ff",
    cell=(40, 40, 80),
    periodic=(True, True, False),
    region=region,
)
```

対応する構造形式は XYZ、MOL2、LAMMPS data です。`Slab` は既存構造を読み込む
オブジェクトで、結晶面自体を生成する API ではありません。スラブは counts・比率の
対象外ですが、`total_atoms` にはスラブ原子も含まれます。

## System

主な属性とメソッド:

```python
system.n_atoms
system.n_chains
system.n_molecules
system.n_packed_molecules
system.positions
system.atom_types
system.masses
system.charges
system.bonds
system.angles
system.dihedrals
system.impropers
system.identity
system.forces
system.box
system.composition
system.requested
system.report
system.write_lammps("out", prefix="system")
system.write_xyz("system.xyz")
```

配列は独立した Python リストです。インデックスは0始まりで、該当しない identity は
`-1` です。リストを書き換えてもネイティブ `System` 自体は変更されません。

## 共通オプション

- `cell=(Lx, Ly, Lz)` または `density=...` で体積を指定します。
- `cell=(Lx, Ly, None), density=...` は残る辺を密度から解きます。
- 長さは Å、密度は g/cm³、温度は K です。
- `periodic=(True, True, True)`、`seed=1234`、`sequence_seed=7` が既定です。
- `output_dir=None` なら中間出力を作らず、所有権を持つ `System` を返します。
- `forces=True` で解析的な力を `system.forces` に保持します。

## Force field の解決

`forcefield=` には `.ff` のパス、`register_forcefield()` で登録した名前、または
既知の別名を指定します。公開配布の `pcff` / `pcff_iff` はリポジトリの
`examples/forcefields/` を探索します。別の場所から実行する場合は明示パス、
`POLYSE_DATA_DIR`、`POLYSE_FORCEFIELD_DIR` のいずれかを使用してください。

