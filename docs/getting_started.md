# Getting started

## 1. インストール

Ubuntu 24.04 x86_64 と CPython 3.13 を使用します。リポジトリのルートで:

```bash
python3 -m pip install dist/polypaves-0.3.0-cp313-cp313-linux_x86_64.whl
python3 - <<'PY'
import polypaves
print(polypaves.__version__)
PY
```

仮想環境を使う場合も、その仮想環境の Python が CPython 3.13 である必要があります。
`not a supported wheel on this platform` は Python ABI、OS、CPU のいずれかが不一致です。

## 2. 最初の実行

```bash
python examples/01_pcff_polystyrene/pmma_dp10/build.py
```

`examples/01_pcff_polystyrene/pmma_dp10/output/` に LAMMPS data、style、identity
が作られます。生成された初期構造は平衡化済み材料ではありません。必要な緩和・
平衡化と妥当性確認はシミュレーション条件に合わせて行ってください。

## 3. Python から確認

```python
import polypaves

system = polypaves.build(
    monomer="*CC*",
    dp=10,
    chains=2,
    forcefield="examples/forcefields/opls_alkane.ff",
    density=0.785,
)

print(system)
print(system.box)
print(system.positions[:2])
print(system.composition)
```

入力不整合は `ValueError`、力場や構造ファイルの欠落は `FileNotFoundError`、
型付け・パラメータ不足・配置失敗などは `polypaves.BuildError` です。

## 4. 外部データの検証

```bash
python3 scripts/external_data.py verify pcff
python3 scripts/external_data.py verify iff
```

詳細は [external_data.md](external_data.md) を参照してください。
