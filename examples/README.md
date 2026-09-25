# Python API examples

隣接する開発リポジトリ `projects/paves/examples` の検証入力を、公開 Python API で
同等に表現できる範囲で移植しています。すべて `import polypaves` を使い、ネイティブ
CLI や `.paves` 入力を経由しません。各スクリプトは場所を自動判定するため、
リポジトリ内のどの作業ディレクトリからでも実行できます。

| 分類 | スクリプト | 内容 |
|---|---|---|
| PCFF | [PMMA DP10](01_pcff_polystyrene/pmma_dp10/build.py) | エステルを含む PMMA 1鎖 |
| PCFF | [PS DP10](01_pcff_polystyrene/ps_dp10/build.py) | 芳香族側鎖を持つ PS 2鎖 |
| コポリマー | [明示配列](03_copolymers/explicit_sequence/build.py) | `Copolymer(sequence=...)` |
| 混合系 | [4成分 blend](05_mixtures/polymer_blend/build.py) | 溶媒、2種のポリマー、コポリマー |
| 混合系 | [PCFF solvent blend](05_mixtures/polymer_blend_in_solvent/build.py) | 質量比と目標原子数 |
| 粗視化 | [CG copolymer](07_coarse_grained/cg_copolymer/build.py) | ランダム exact 配列と junction typing |
| 末端 | [非対称末端](10_sequences_and_ends/asymmetric_end_groups/build.py) | head/tail に異なる末端基 |
| PCFF-IFF | [Kapton DP3](12_pcff_iff/kapton_dp3/build.py) | 環構造を含む Kapton |
| 系サイズ | [原子数＋box](13_packing_size/atom_count_explicit_box/build.py) | 単成分の原子数目標 |
| 系サイズ | [blend 原子数](13_packing_size/blend_atom_count/build.py) | mol比と原子数目標 |
| 系サイズ | [鎖数＋box](13_packing_size/chain_count_explicit_box/build.py) | 明示セルと鎖数 |
| 系サイズ | [溶媒＋ポリマー](13_packing_size/solvent_and_polymer_atom_count/build.py) | mol比、原子数目標、明示セル |

例:

```bash
python examples/03_copolymers/explicit_sequence/build.py
```

元の常時検証17例から、Python API で意味を保てる12例を移植しています。ビーズ専用
構文、log-normal 鎖長分布、構築後の座標 perturb は、現在の公開 Python API に対応する
型・引数がないため収録していません。`08_configuration_file` は設定ファイルと CLI
override 自体を確認する例であり、Python API ではその経路を使用しないため対象外です。
異なる意味へ置き換えた見かけだけの変換は行っていません。

`examples/forcefields/` の OPLS/CG 設定は例の一部です。PCFF/IFF 設定は
`external/` のパラメータと typing template を参照します。生成物は各例の
`output/` に保存され、Git には追加されません。
