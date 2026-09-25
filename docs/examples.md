# Examples

公開例はすべて Python API を直接呼び出します。

```text
build.py
  └─ import polypaves
       └─ polypaves.build(...) / polypaves.pack(...)
            └─ native C++ engine
                 └─ System
```

`.paves` 設定ファイルを解釈する CLI をサブプロセス起動する構成ではありません。
Python が所有する `System` が返り、出力前に座標、型、電荷、結合、組成、box、
構築 report を検査できます。

全例の一覧、元の開発例との対応、API の機能差または CLI 固有の目的により移植していない例は
[examples/README.md](../examples/README.md) に記載しています。

単一成分では `polypaves.build()`、複数成分・コポリマー・スラブでは
`polypaves.pack()` を使います。`counts` と組成比＋目標サイズは排他的です。

```python
import polypaves

ps = polypaves.Polymer("*C(c1ccccc1)C*", dp=20, name="ps")
pmma = polypaves.Polymer("*CC(C)(C(=O)OC)*", dp=20, name="pmma")

system = polypaves.pack(
    [ps, pmma],
    total_molecules=20,
    mole_fractions={"ps": 0.5, "pmma": 0.5},
    forcefield="pcff",
    density=1.0,
)
print(system.composition)
```
