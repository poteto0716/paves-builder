# network_propertiesで計算する物性

単分子プロジェクト（`monomer`）と架橋プロジェクト（`crosslink`）で、1レプリカごとに計算する物性の一覧です。
値は各レプリカのJSONと、プロジェクトの`property_summary.csv`／`replica_properties.csv`／`report.md`に入ります。
架橋密度ごとの比較は単分子プロジェクトの`crosslink_series.csv`／`.png`と`report.md`の「Crosslink series」に入ります。
レプリカ間の平均・SD・SEMは、チェックを通らなかった値も含む`mean`と、通った値だけの`checked_mean`の両方を出します。

## 計算の流れ

```
単分子（monomer）
  packing（0.5 g/cm³）→ T・200 MPa加圧NPT（50 ps）→ T・1 atm NPT（1 ns）        … 開始点なし・反応なし
  → 熱履歴：300→800 K加熱、800 K緩和、800→100 K冷却（40 K/ns）            … Tg・CTE
  → 100→300 K、300 K緩和 → 300 K NPT（10 ns）、300 K NVT（10 ns）          … 300 K熱力学量・環の数密度
  → 300 K→T、T緩和（1 ns）→ T NPT（2 ns）→ T NVT（2 ns）                   … Tの熱力学量・環の統計
     └ TのNPTの最後の構造 = 架橋の出発構造

架橋（crosslink、密度ごと）：レプリカNは単分子レプリカNの出発構造から
  開始点作成 → opt → T・NVT 0.1 fs（1 ps）→ T・NPT（100 ps）
  → Tで成長cycle（目標反応率99.5%）→ 最終緩和 →（任意）ゾル削除             … 網目構築の量
  → 熱履歴（Tg・CTE）→ 300 K NPT・NVT                                       … 300 K物性
  → 網目グラフ解析（300 K構造）、溶媒（指定時）
  → 弾性（任意：mechanicalコマンド）
```

Tは既定363.15 K（90 ℃）。温度・時間は既定値です。各量の詳しい定義は`README.md`と`../crosslinking/README.md`にあります。

---

## 1. 単分子の構築（Monomer build）

出典: `replica_XX/build/done.json`

| 物性 | キー | 意味 |
|---|---|---|
| 分子数・原子数 | `build.molecules`／`build.atoms` | packingした分子数（総原子数の目標既定20000を超える最小の分子数） |
| 環の数 | `build.reactive_sites` | 開始点ルールの`target_smarts`に一致する環の数（エポキシ・オキセタンなど） |

## 2. 単分子のTでの物性（<T> K monomer liquid）

出典: `replica_XX/<T>K/result.json`（環の統計の詳細は`<T>K/03_production/reactive_rings/`）。
熱履歴と300 Kの段階の後、300 K NPTの続きからTへ昇温し、1 ns緩和した後のNPT 2 ns（1 psごと）とNVT 2 nsから求めます。
表の「Epoxide」「Oxetane」はchemistryに合わせた環の名前です。

| 物性 | キー | 求め方・意味 | 注意 |
|---|---|---|---|
| 密度・Cp・Cv・圧縮率・体積弾性率・CTE・音速（T） | `at_T.npt.*` | 「300 Kの熱力学量」と同じ揺らぎの式をTのNPTに適用 | |
| 環の数密度（T） | `at_T.rings.number_density_nm3`／`_mmol_cm3` | 環の数／TのNPT平均体積。1 nm⁻³ = 1/0.6022 mmol/cm³ | |
| 配位数（4 Å、6 Å） | `at_T.rings.coordination_<r>A` | 1つの環あたりの接触相手の平均数。一方の環のOと他方の環の攻撃されるCH2炭素の距離がカットオフ以内なら接触 | リレー反応1段の幾何。カットオフは反応の開始距離と上限（`ring_contact_cutoffs_A`）。同一分子内も含む |
| 接触グラフの最大クラスター分率（4 Å、6 Å） | `at_T.rings.largest_cluster_<r>A` | 接触グラフの最大連結成分に入る環の数／全環数 | 反応の連鎖がどこまで届きうるかの目安 |
| 環のMSD | `at_T.rings.msd_final_A2` | 環の重原子重心の平均二乗変位（観測終了時のラグ）、Å²。周期境界をまたいで展開し、系の重心ドリフトを除き、全時間原点で平均 | |
| 環の拡散係数 | `at_T.rings.diffusion_coefficient_m2_s` | ラグ10〜50%の範囲でMSD/6t、m²/s | MSD指数αが0.8〜1.2（拡散域）のときだけチェック通過 |
| MSD指数α | `at_T.rings.msd_exponent_alpha` | 同じ範囲のlog(MSD)–log(t)の傾き | 1なら拡散域、1未満はケージ内の運動が残っている |
| 直接Cv・誘電応答（T） | `at_T.nvt.*` | 「300 K NVT」と同じ式をTのNVTに適用 | |

環の定義はエポキシ・オキセタン共通です（環のO、Oに結合したH2を持つ環炭素 = 攻撃されるCH2、環の重原子）。
SMARTSが同じ環に複数回一致しても（オキセタンは環のCH2のどちらも1番目に当たる）1つと数えます。

## 3. 網目構築（Network construction、架橋のみ）

出典: `replica_XX/network_origin.json`（crosslinkingの`metrics.json`から転記）

| 物性 | キー | 求め方・意味 | 注意 |
|---|---|---|---|
| 架橋反応率 | `crosslinking.conversion` | 開始点作成後に残っていた環のうち、成長反応で開いた割合 | 目標（既定0.995）に届かないレプリカはチェック不合格 |
| 開始点を含む環の消費率 | `crosslinking.epoxide_consumption_including_activation` | （開始点数 + 成長反応数）／初期の環の数 | 系全体（ゾル込み） |
| 開始点数・開始点の割合 | `crosslinking.activation_count`／`actual_activation_fraction` | 実際に作った開始点の数と、環の数に対する割合（割合指定は環の数×割合の切り捨て） | |
| 体積あたりの開始点数 | `crosslinking.activation_density_nm3`／`_mmol_cm3` | 開始点数／架橋開始直前（開始点作成後のNPT緩和の後）のセル体積 | |
| 活性末端あたりの連鎖長（数平均・重量平均） | `crosslinking.kinetic_chains.*_average_epoxide_units` | 反応履歴で活性（OH*）の移動をたどり、各開始点が何個の環を開いたか（開始点の環 + 成長反応数）。数平均と重量平均 | 架橋点間の鎖長（下の「網目トポロジー」）とは別物 |
| 連鎖長の分散度・最長 | `crosslinking.kinetic_chains.dispersity`／`max_epoxide_units` | 重量平均／数平均、最長の連鎖 | 分布は`chain_lengths.csv`／`.png` |
| ゲルの環の反応率・ゾル分率 | `crosslinking.gel_epoxide_conversion`／`sol_mass_fraction` | ゾル削除後のゲルについて、反応済の環／ゲルのモノマーがもともと持っていた環（開始点で開いた環を含む）。削除したゾルの質量分率 | `--remove-sol-fraction`のときだけ |

ゾル削除（`--remove-sol-fraction`）では、周期次元3の分子だけを残します。周期次元は、結合をたどって同じ原子の別の周期像に行き着いたときの
セルオフセット（周期ベクトル）の階数で、0は有限分子、1は紐、2はシート、3は3次元網目です。斜めに巻き付いた紐を3方向貫通と誤認しません。
ゾル分率と削除数は`network_origin.json`に入ります。

## 4. Tgと熱膨張（Cooling）

出典: `dry_result.json`の`cooling`。800→100 Kを40 K/nsで冷却し、各温度plateauの後半の密度を使います。

| 物性 | キー | 求め方・意味 | 注意 |
|---|---|---|---|
| Tg | `cooling.tg_K` | 10 K幅で平均した密度–温度を、2本の直線を滑らかにつないだ双曲線でフィットした折れ曲がり点 | 冷却速度に依存する見かけ値（実験より高く出る） |
| Tgでの密度 | `cooling.density_at_tg_g_cm3` | フィット曲線のTgでの密度 | |
| 線CTE（glass／rubber） | `cooling.alpha_linear_glass_1_K`／`_rubber_1_K` | フィット曲線のTg−50 K／Tg+50 Kでの解析微分dρ/dTと密度ρ(T)から αL = −(dρ/dT)/(3ρ(T)) | 遷移幅が広いとglassにrubberの傾きが混ざる（`cte_mixing_fraction`）。Tg±50 Kが冷却範囲外ならチェック不合格 |
| 体積CTE（glass／rubber） | `cooling.alpha_volume_*_1_K` | 3 αL | 等方性を仮定 |

チェック（`fit_checks_passed`）の条件は、R² > 0.95、遷移幅 < 温度範囲の1/4、直線よりAICが6以上改善、
αL(glass) > 0、Tgが範囲の端10%に張り付いていない、Tg±50 Kが範囲内、高温での沸騰（密度が開始時の0.7倍未満）がない、の7つです。
`report`／`collect`の時に、古い版のフィットは`cooling.csv`から再フィットされます。

## 5. 300 Kの熱力学量（300 K thermodynamics）

出典: `dry_result.json`の`npt`。熱履歴後の300 K・1 atm NPT観測（10 ns、1 psごと）の揺らぎから求めます。

| 物性 | キー | 求め方 | 注意 |
|---|---|---|---|
| 密度 | `npt.density_g_cm3` | 質量／平均体積、g/cm³ | |
| 定圧比熱 Cp | `npt.Cp_J_kg_K` | Var(H)/(kB T² m) | 古典・拘束付きSageモデル。量子補正なし（実測より大きく出やすい） |
| 定積比熱 Cv（NPT恒等式） | `npt.Cv_J_kg_K` | Cp − T αV²/(ρ κT) | |
| 等温圧縮率 κT／等温体積弾性率 KT | `npt.compressibility_isothermal_1_Pa`／`bulk_modulus_isothermal_GPa` | Var(V)/(kB T⟨V⟩)、KT = 1/κT | |
| 断熱圧縮率 κS／断熱体積弾性率 KS | `npt.compressibility_isentropic_1_Pa`／`bulk_modulus_isentropic_GPa` | κS = κT Cv/Cp、KS = 1/κS | |
| 体積CTE αV（300 K） | `npt.alpha_volume_1_K` | Cov(V, U+PV)/(kB T²⟨V⟩) | 揺らぎから求める300 Kの値。冷却曲線のCTEとは別の求め方 |
| 線CTE αL（300 K） | `npt.alpha_linear_1_K` | αV/3 | |
| bulk音速 | `npt.bulk_sound_speed_m_s` | √(KS/ρ) | 固体の縦波音速とは異なる |
| 環の数密度（300 K、単分子のみ） | `rings.number_density_300K_nm3`／`_mmol_cm3` | 環の数／300 K NPT平均体積 | 架橋膜の未反応環は「網目トポロジー」のUnreacted epoxide／oxetane |

5ブロックの誤差、有効サンプル数、前後半のdriftでサンプリングのチェックをします。

## 6. 300 Kの構造（300 K structure）

| 物性 | キー | 求め方・意味 | 注意 |
|---|---|---|---|
| 残膜率（架橋のみ） | `structure.film_retention_volume_ratio` | V(300 K、架橋・熱履歴後のNPT平均)／V(300 K、同じ番号の単分子レプリカのNPT平均) | 両方300 Kなので熱膨張は入らない。ゾル削除時は分子がゲルだけの体積。現像による膜減りの目安で、膨潤・乾燥収縮や深さ方向の露光むらは含まない |
| 点プローブ空隙率 | `structure.point_probe_void_fraction` | NPT観測後のスナップショットにランダム点（10000点）を置き、どの重原子のvan der Waals半径（C 1.7 Å、O 1.52 Åなど）の内側にも入らない点の割合。Hは無視 | 時間平均でも、溶媒がアクセスできる体積でもない |

## 7. 網目トポロジー（Network topology, 300 K体積、架橋のみ）

出典: `dry_result.json`の`network`。300 K NPT観測後の構造の化学グラフから求めます（MD不要）。
H原子を除いた重原子グラフで、3〜8員環の環系（縮環を含む）は1つの節点にまとめます。

| 物性 | キー | 意味 |
|---|---|---|
| 有効鎖密度 ν_e | `network.strand_density_nu_e_mmol_cm3` | 架橋点どうしを結ぶ有効な鎖の数密度 |
| 有効鎖の本数 | `network.effective_strands` | |
| Mc（セル全質量／鎖本数） | `network.Mc_g_mol` | ρN_A/ν_e。ゾル・ダングリング・ループの質量も含む見かけの値 |
| Mc（有効網目の質量／鎖本数） | `network.Mc_active_network_g_mol` | 架橋点の原子を含む |
| 鎖の平均分子量 | `network.strand_molecular_weight_g_mol.mean` | 架橋点の間の原子（H込み）の分子量。架橋点の原子は含まない |
| 鎖の平均長さ | `network.strand_heavy_atoms.mean` | 架橋点の間の重原子数（途中の環系は全原子を数える） |
| 鎖の主鎖結合数 | `network.strand_backbone_bonds.mean` | 架橋点の原子間の最短経路の結合数 |
| 架橋点の数・密度・平均官能度 | `network.junctions`／`junction_density_mmol_cm3`／`mean_junction_functionality` | 有効な結合を3本以上持つ節点 |
| cycle rank密度 ξ/V | `network.cycle_rank_density_mmol_cm3` | ξ = 鎖の数 − 架橋点の数 |
| ダングリング質量分率 | `network.dangling_mass_fraction` | 1本の結合だけでつながり、その先に無限につながった部分がない部分 |
| ループ質量分率 | `network.loop_mass_fraction` | 同じ架橋点に戻る弾性的に無効なループ |
| ゾル質量分率 | `network.sol_mass_fraction` | 周期次元0（有限）の分子の質量分率。周期次元1以上をゲルとみなす基準 |
| 3次元ゲル質量分率 | `network.three_dimensional_gel_mass_fraction` | 周期次元3の分子の質量分率（ゾル削除と同じ厳しい基準）。次元ごとの内訳は`mass_fraction_by_periodic_dimension` |
| ゴム状弾性率の目安（affine／phantom） | `network.rubbery_shear_modulus_estimates_MPa.*` | ν_e kT／ξkT/V（300 K）。からみ合いは含まない。Tgより上のゴム状態の目安で、300 Kのガラスの弾性率とは別物 |
| 極性基の密度 | `network.polar_groups_mmol_cm3.*` | 脂肪族OH、フェノール性OH、エーテルO、未反応エポキシ環・オキセタン環、芳香環（mmol/cm³）。OH*の代理表現では活性末端もOHとして数える |

鎖の長さ・分子量の全分布は`network_strands.csv`／`.png`に出ます。

## 8. 300 K NVT（300 K NVT）

出典: `dry_result.json`の`nvt`。300 K NPT観測後の体積で10 ns NVT。

| 物性 | キー | 求め方 | 注意 |
|---|---|---|---|
| 直接Cv | `nvt.Cv_direct_J_kg_K` | Var(E)/(kB T² m) | NPT恒等式のCvと比べる検証用 |
| 誘電応答（固定電荷）、等方平均とテンソル | `nvt.dielectric_isotropic_fixed_charge`／`nvt.dielectric_xx`… | 全電荷変位（双極子）の揺らぎから、PME導体境界の揺らぎ公式で求める | 電子分極を含まない（実測の高周波分を欠く）。ガラスの配向緩和が観測時間より長いと収束しないので、自動では収束済みと判定しない |

## 9. 弾性（Elasticity、任意）

既定では計算しません。架橋膜ができた後に`mechanical <プロジェクト>`（または`submit --mechanical`）で実行します。
出典: `replica_XX/mechanical/result.json`。300 K構造から6成分×正負1%のひずみへ分岐し、各200 ps NVT緩和・1 ns観測。
応力は全ポテンシャルの微小セル変形差分（PME込み）と運動項から求めます。X–H結合は拘束せず調和結合で扱い、0.25 fs刻みです。

| 物性 | キー | 意味 |
|---|---|---|
| 体積弾性率（Voigt／Reuss／Hill） | `mechanical.bulk_*_GPa` | 6×6剛性行列からの等方平均 |
| ずり弾性率（Voigt／Reuss／Hill） | `mechanical.shear_*_GPa` | 同上 |
| ヤング率（Hill） | `mechanical.young_Hill_GPa` | |
| ポアソン比（Hill） | `mechanical.poisson_Hill` | |
| 剛性行列 C | `mechanical.C_<ij>_<kl>_GPa` | 6×6（CSVのみ） |

行列に非正の固有値があれば等方弾性率は確定しません。ひずみ幅・応力の誤差・セルサイズ依存の確認が必要なため、
自動の収束フラグは常にfalseです（比較用の推定値）。

## 10. 溶媒（Solvent、架橋で指定時のみ）

出典: `replica_XX/solvents/<name>/result.json`。300 Kの網目に対して溶媒ごとに計算します。

| 物性 | キー | 求め方・意味 | 注意 |
|---|---|---|---|
| 挿入自由エネルギー | `solvents.<name>.coupling_free_energy_kJ_mol` | 乾燥樹脂へ1分子をsoft-core LJと電荷の22窓で段階的に挿入し、MBARで求める | 孤立分子（気相）基準の過剰自由エネルギー。液体からの移行自由エネルギーではない |
| 相対化学ポテンシャル | `solvents.<name>.chemical_potential_relative_kJ_mol` | 上に有限セルの濃度項（参照1 mol/L）を加えた量 | |
| 溶媒の自己拡散係数 | `solvents.<name>.diffusion_coefficient_m2_s` | 指定数の溶媒をpackし、10 ns NPTのMSD/6t | 指定数は平衡吸収量ではない。拡散域に届かなければnull |

## 注意（全体）

- 力場は固定電荷のOpenFF Sage、古典MDです。比熱の量子補正、電子分極、反応の活性化障壁は含みません。
- 架橋は距離基準の反応（開始4 Å→最大6 Å）で、反応速度論は再現しません。連鎖長や反応率は構築手順に依存します。
- 2万原子・数nsのMDなので、Tg（冷却速度）、弾性（ひずみ・観測時間）、誘電応答（配向緩和）は特に有限時間・有限サイズの影響を受けます。
- スモーク実行（`--smoke`）の値に物理的な意味はありません（チェックは常に不合格）。
