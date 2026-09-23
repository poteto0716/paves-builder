# External data

このリポジトリの `external/` には examples を再現するための固定版 PCFF/IFF データと
silica 座標を収録しています。これらは polyse の所有物ではなく、polyse の
[LICENSE](../LICENSE) は適用されません。利用・再配布条件は各権利者の条件を確認して
ください。以下の手順は、収録コピーを公式配布物から checksum 検証付きで復元する
ためのものです。EMC executable は polyse に使用・リンクしていません。

## PCFF polymer examples

Obtain the EMC distribution from the
[author's download page](https://montecarlo.sourceforge.net/emc/Welcome.html) or
[official file listing](https://sourceforge.net/projects/montecarlo/files/).
The validated files in the official `emc_linux_x86_64_v9.4.4_20260701.tgz`
archive are under `v9.4.4/field/cff/pcff/` (some older/local installations
expose them as `field/pcff/`):

- `pcff.frc` (header includes version 5.1, 1 December 2020)
- `pcff_templates.dat`

After unpacking the distribution, pass the directory containing these two files:

```bash
python3 scripts/external_data.py pcff /path/to/v9.4.4/field/cff/pcff
```

The official archive was downloaded during package validation and both database
files matched the pinned hashes. This copies only the two checksum-verified files into `external/pcff/`.
The absolute path above is your one-time installation input; example configurations
use only relative paths within this repository.

The inspected EMC installation includes GPLv3 license text, while its SourceForge
project metadata labels the project MIT. Neither observation by itself establishes
rights to redistribute the separately sourced historical PCFF database.収録されていること自体
からライセンスを推定したり、polyse が権利を付与したりするものではありません。

## IFF silica/polymer interface

Use the [official INTERFACE MD page](https://bionanostructures.com/interface-md/)
and its release 1.5 link, also referenced by the
[authors' repository](https://github.com/hendrikheinz/INTERFACE-force-field-and-surface-models).
The exact archive used for validation is
[interface_ff_1_5.zip](https://bionanostructures.com/wp-content/uploads/2016/02/interface_ff_1_5.zip).

After obtaining that archive:

```bash
python3 scripts/external_data.py iff /path/to/interface_ff_1_5.zip
```

Archive SHA-256:

```text
42b26c36ee91254eb83b2e87f8ac4dd8258dd59b4d7aa966c5fd8123410a2086
```

Only these four archive members are installed under `external/iff/`:

| Archive member below `INTERFACE_FF_1_5/` | Installed filename |
|---|---|
| `FORCE_FIELDS/pcff_interface_v1_5.frc` | `pcff_interface_v1_5.frc` |
| `FORCE_FIELDS/pcff_interface_v1_5_templates.dat` | `pcff_interface_v1_5_templates.dat` |
| `MODEL_DATABASE/SILICA/silica_Q3_amorph_4_7OH_0pct_ion.car` | `silica.car` |
| `MODEL_DATABASE/SILICA/silica_Q3_amorph_4_7OH_0pct_ion.mdf` | `silica.mdf` |

必要に応じて CAR/MDF から生成した MOL2 や LAMMPS 出力は Git 管理対象外です。
元の PCFF/IFF と CAR/MDF の固定コピーだけを `external/` に収録しています。

## Verification and citations

```bash
python3 scripts/external_data.py verify pcff
python3 scripts/external_data.py verify iff
```

Individual hashes are in `scripts/external_checksums.json`. A mismatch is refused;
a newer file requires a separate validation, not disabling the checksum.

For IFF cite H. Heinz, T.-J. Lin, R. K. Mishra, F. S. Emami,
*Langmuir* 2013, 29, 1754, [DOI](https://doi.org/10.1021/la3038846).
For the silica models cite F. S. Emami et al., *Chem. Mater.* 2014, 26, 2647,
[DOI](https://doi.org/10.1021/cm500365c).
Consult the PCFF database's reference records for the original parameter sources.

第三者データは polyse wheel 自体には格納していません。`external/` のファイルを
別の成果物へ再梱包する場合も、各権利者の条件を別途確認してください。
