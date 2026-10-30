"""SMILES／総原子数からpure epoxyをpackし、初期NPT緩和後に反応MDを行う。

packing後の加圧NPT: 300 K、200 MPa、1 fs × 50,000 step = 50 ps（void除去）。
初期NPT: 300 K、1 atm、1 fs × 1,000,000 step = 1 ns。

各cycle: 反応 → opt → NVT (300 K, 0.1 fs, 10000 step)
                    → NPT (300 K, 1 atm, 1 fs, 1000 step)。
反応後の緩和は--relaxation-stages（JSON文字列かJSON/YAMLファイル）で任意の段に置き換えられる。
OpenMMはCUDA/mixedが既定。CPUでは --platform CPU を指定する。
OH*は中性graphにcationic_activityを付けた既存のproxy表現。
OVITO軌跡とcycle別構造は既定では保存しない。--save-trajectoryで有効にする。
保存先はカレントのresults/<project>。既存projectの上書きは --force が必要。
activation.yaml／reactions.yamlはカレントディレクトリから読み込む。
反応距離は5 Åから開始。反応0件が2 cycle続くと1.2倍に拡大し、最大7 Å（--max-reaction-cutoff）。
反応後の割り当ては局所Sageフラグメントの電荷テンプレートを再利用する。
"""
import argparse
import copy
from contextlib import nullcontext
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import re
import shutil
from types import SimpleNamespace
import polypaves as p
from polypaves.structure import write_structure_files
from polypaves.reaction import (
    ActivationInitializer, ReactiveMDRunner, MDStage, AdaptiveReactionCutoff, LocalSageForceField,
    OpenMMBridge, OpenMMStateSnapshot,
    OVITOTrajectoryWriter,
    load_activation_definitions, load_reaction_definitions,
)

DEFAULT_SMILES = 'CCC(COCC1CO1)(COCC2CO2)COCC3CO3'
# TMP triglycidyl etherはH込み47原子/分子。開始点1個につきHを2個追加する。
DEFAULT_TOTAL_ATOMS = 20_000
DEFAULT_ATOM_COUNT_MODE = 'above'
DEFAULT_CYCLES = 5000
DEFAULT_TARGET_CONVERSION = .995
STAGE_KEYS = {'ensemble', 'steps', 'timestep', 'temperature', 'pressure_atm'}


def periodic_winding(graph, xyz, cell):
    """Axes along which one molecule wraps onto its own periodic image.

    Walks the bonds with minimum-image displacements; reaching an atom again at
    a different cell offset means the molecule is infinite along that axis.
    xyz and cell are in angstrom; cell rows are the box vectors.
    """
    import numpy as np
    inverse = np.linalg.inv(cell)
    neighbors = [[] for _ in graph['atoms']]
    for i, j, _ in graph['bonds']:
        shift = -np.rint((xyz[j] - xyz[i]) @ inverse).astype(int)
        neighbors[i].append((j, shift))
        neighbors[j].append((i, -shift))
    lifts, winding = {}, np.zeros(3, dtype=bool)
    for root in range(len(neighbors)):
        if root in lifts:
            continue
        lifts[root], pending = np.zeros(3, dtype=int), [root]
        while pending:
            i = pending.pop()
            for j, shift in neighbors[i]:
                expected = lifts[i] + shift
                if j not in lifts:
                    lifts[j] = expected
                    pending.append(j)
                else:
                    winding |= expected != lifts[j]
    return winding.tolist()


def without_sol(structure, origin):
    """Copy of structure keeping only molecules that percolate through the cell.

    Development dissolves every finite molecule: unreacted monomers and
    oligomers not attached to the infinite network (the sol fraction). A
    molecule is kept when it wraps onto its own periodic image along at least
    one axis. origin maps stable_atom_id -> index of the initial (activated)
    monomer, to report removed monomers and oligomers separately.
    Returns (structure or None when nothing percolates, summary dict).
    """
    import numpy as np
    xyz, cell = np.asarray(structure._positions), np.asarray(structure.box_vectors)
    offsets = np.cumsum([0] + [len(g['atoms']) for g in structure.graphs])
    masses = [sum(a['mass'] for a in g['atoms']) for g in structure.graphs]
    keep, axes, monomers, oligomers = [], [], 0, 0
    for i, graph in enumerate(structure.graphs):
        winding = periodic_winding(graph, xyz[offsets[i]:offsets[i + 1]], cell)
        if any(winding):
            keep.append(i)
            axes.append(winding)
        elif len({origin[a['stable_atom_id']] for a in graph['atoms'] if a['stable_atom_id'] in origin}) == 1:
            monomers += 1
        else:
            oligomers += 1
    kept_atoms = int(sum(offsets[i + 1] - offsets[i] for i in keep))
    summary = dict(sol_monomers_removed=monomers, sol_oligomers_removed=oligomers,
                   sol_atoms_removed=structure.n_atoms - kept_atoms,
                   sol_mass_fraction=1 - sum(masses[i] for i in keep) / sum(masses),
                   network_molecules_kept=len(keep), network_winding_axes=axes, final_total_atoms=kept_atoms)
    if not keep:
        return None, summary
    trimmed = copy.copy(structure)
    atoms = np.concatenate([np.arange(offsets[i], offsets[i + 1]) for i in keep])
    trimmed.graphs = [structure.graphs[i] for i in keep]
    trimmed._components = [structure._components[i] for i in keep]
    trimmed._positions = xyz[atoms].copy()
    if hasattr(structure, '_original_identity'):
        trimmed._original_identity = [structure._original_identity[a] for a in atoms]
    trimmed.report = {**structure.report, 'atoms': trimmed.n_atoms, 'molecules': len(trimmed.graphs)}
    return trimmed, summary


def load_relaxation_stages(value):
    """JSON文字列、またはJSON/YAMLファイルのパスから反応後緩和の段リストを読む。

    例: [{"ensemble": "NVT", "steps": 5000, "timestep": 0.1},
         {"ensemble": "NPT", "steps": 20000, "timestep": 1.0, "pressure_atm": 1}]
    ファイルはリストそのもの、または{stages: [...]}。
    """
    if value is None or isinstance(value, list):
        return value
    import yaml
    path = Path(value)
    data = yaml.safe_load(path.read_text() if path.is_file() else value)
    if isinstance(data, dict) and set(data) == {'stages'}:
        data = data['stages']
    return data


def cycle_stages(spec, *, temperature, pressure_atm):
    """段リストをMDStageへ変換する。timestep: fs、temperature: K（省略時は全体の温度）、
    pressure_atm: NPTのみ（省略時は全体のpressure_atm）。"""
    if not isinstance(spec, list) or not spec:
        raise ValueError('relaxation_stages must be a non-empty list of stages')
    stages = []
    for index, entry in enumerate(spec, 1):
        if not isinstance(entry, dict) or not {'ensemble', 'steps', 'timestep'} <= set(entry) \
                or set(entry) - STAGE_KEYS:
            raise ValueError(f'relaxation stage {index} needs ensemble, steps, timestep and may add '
                             f'temperature, pressure_atm: {entry!r}')
        ensemble = str(entry['ensemble']).upper()
        if type(entry['steps']) is not int or entry['steps'] < 1:
            raise ValueError(f'relaxation stage {index}: steps must be a positive integer')
        if ensemble == 'NVT' and 'pressure_atm' in entry:
            raise ValueError(f'relaxation stage {index}: NVT must not set pressure_atm')
        pressure = entry.get('pressure_atm', pressure_atm) * 1.01325 if ensemble == 'NPT' else None
        stages.append(MDStage(ensemble, steps=entry['steps'], timestep=entry['timestep'],
                              temperature=entry.get('temperature', temperature), pressure=pressure))
    return tuple(stages)


def activation_yaml_with_selection(text, *, count=None, fraction=None):
    """Return activation.yaml text whose first rule uses the CLI count/fraction.

    Only the selection/count/fraction lines of the first rule (the one this
    script reads) are rewritten, so comments and other keys stay as written.
    With neither override the text is returned unchanged.
    """
    if count is None and fraction is None:
        return text
    lines = text.splitlines(keepends=True)
    starts = [i for i, line in enumerate(lines) if re.match(r'\s*- ', line)]
    if not starts:
        raise ValueError('activation.yaml has no activation rule')
    first = starts[0]
    item = len(lines[first]) - len(lines[first].lstrip())
    end = next((i for i in starts[1:] if len(lines[i]) - len(lines[i].lstrip()) <= item), len(lines))
    key = re.compile(r'(\s*)(selection|count|fraction):')
    found = [(i, key.match(lines[i])) for i in range(first, end) if key.match(lines[i])]
    selections = [(i, m) for i, m in found if m.group(2) == 'selection']
    if len(selections) != 1:
        raise ValueError('the first activation rule must have exactly one selection line')
    index, match = selections[0]
    indent, newline = match.group(1), '\n'
    mode, value = ('count', count) if count is not None else ('fraction', fraction)
    replacement = [f'{indent}selection: {mode}{newline}', f'{indent}{mode}: {value!r}{newline}']
    drop = {i for i, m in found if m.group(2) in ('count', 'fraction') and m.group(1) == indent}
    body = [line for i, line in enumerate(lines[first:end], first) if i not in drop and i != index]
    body.insert(index - first - sum(i < index for i in drop), ''.join(replacement))
    return ''.join(lines[:first] + body + lines[end:])


def plan_size(total_atoms, atoms_per_molecule, sites_per_molecule, *, added_per_site,
              count=None, fraction=None, mode=DEFAULT_ATOM_COUNT_MODE):
    """Resolve whole molecules including activation atoms in the final budget."""
    if mode not in ('above', 'nearest', 'floor', 'exact'):
        raise ValueError('atom_count_mode must be above, nearest, floor or exact')
    if (count is None) == (fraction is None):
        raise ValueError('Specify one activation count or fraction')
    if count is not None and (type(count) is not int or count < 1):
        raise ValueError('activation_count must be a positive integer')
    if fraction is not None and (isinstance(fraction, bool) or not math.isfinite(fraction)
                                 or not 0 < fraction < 1):
        raise ValueError('activation_fraction must be between 0 and 1, excluding endpoints')

    def size(n):
        activations = count if count is not None else math.floor(fraction * (sites_per_molecule * n))
        return n * atoms_per_molecule + added_per_site * activations, activations

    lower, upper = 0, total_atoms // atoms_per_molecule + 1
    while upper - lower > 1:
        midpoint = (lower + upper) // 2
        if size(midpoint)[0] <= total_atoms:
            lower = midpoint
        else:
            upper = midpoint
    molecules = lower
    if mode == 'above':
        molecules = upper
    elif mode == 'nearest':
        molecules = min((lower, upper), key=lambda n: (abs(size(n)[0] - total_atoms), -n))
    actual_atoms, activations = size(molecules)
    if mode == 'exact' and actual_atoms != total_atoms:
        raise ValueError(f'exact total_atoms is unreachable; neighboring totals are {size(lower)[0]} and {size(upper)[0]}')
    if molecules < 1 or activations < 1 or activations >= sites_per_molecule * molecules:
        raise ValueError('Size must allow at least one start and one remaining epoxide')
    return dict(monomers=molecules, activation_count=activations, actual_total_atoms=actual_atoms,
                packed_atoms=molecules * atoms_per_molecule,
                initial_epoxide_sites=molecules * sites_per_molecule)


def run(project, *, force=False, compression_steps=50_000, compression_pressure_mpa=200.,
        equilibration_steps=1_000_000, smiles=DEFAULT_SMILES, total_atoms=DEFAULT_TOTAL_ATOMS,
        atom_count_mode=DEFAULT_ATOM_COUNT_MODE, density=.5, activation_count=None, activation_fraction=None,
        platform='CUDA', cycles=DEFAULT_CYCLES, temperature=300.,
        equilibration_timestep=1.,
        nvt_timestep=.1, nvt_steps=10000, npt_timestep=1., npt_steps=1000,
        pressure_atm=1., minimization_iterations=2000,
        target_conversion=DEFAULT_TARGET_CONVERSION, max_reactions_per_cycle=None, seed=31415,
        reaction_cutoff=5., max_reaction_cutoff=7., cutoff_growth_factor=1.2,
        cutoff_zero_cycles=2, fixed_reaction_cutoff=False,
        relaxation_stages=None, post_reaction_md_steps=0, stop_after_idle_cycles=200,
        reaction_minimization_radius=12., final_relaxation_steps=1000, final_md_steps=1_000_000,
        pack_relax=50, pack_postrelax=400,
        template_workers=4, reaction_assignment='local_templates', local_radius=5, local_context_radius=8,
        max_fragment_atoms=1200, template_cache=None, require_reaction_templates=False,
        save_trajectory=False, remove_sol_fraction=False):
    """packingから反応MD、最終構造の保存まで実行する。

    project: 保存先はカレントディレクトリ/results/<project>。
      名前は英数字・Unicode文字・_で始まり、以降はこれらと.-を使用できる。
    force: Trueなら既存projectディレクトリの内容を削除して新規実行する。
    activation.yaml／reactions.yamlはカレントから読み、保存先にもコピーする。
    total_atoms: 初期化後の目標総原子数（H込み、既定20000）。above（既定）は
      目標を厳密に超える最小の整数分子数、nearestは最も近い整数分子数、
      floorは目標以下、exactは厳密一致（到達不能ならエラー）。
    density: packing時のg/cm^3（既定0.5）。0.8で詰めた嵩高い芳香族エポキシは、200 MPa・1 fsの加圧NPTで
      半数がNaNになった（0.5では0/4）。加圧後は約1.1 g/cm^3まで縮む。temperature: K。timestep: fs。
    pressure_atm: atm（ライブラリのbarへ換算して渡す）。
    compression_steps: 加圧NPTのstep数（既定50,000）。equilibration_timestep刻みで既定50 ps。
      0なら加圧しない。compression_pressure_mpa: 加圧NPTの圧力 MPa（既定200）。
      packing直後のvoidを潰すため、初期NPTより前に実行する。
    equilibration_steps: 初期NPTのstep数（既定1,000,000）。既定1 fsで1 ns。
    packing → 開始点作成 → opt → 加圧NPT → 初期NPT緩和 → 成長cycle。
    各cycleは反応探索／graph更新／Sage再割当／opt → NVT → NPT。
    cycles: 最大反応サイクル数（既定5000）。target_conversion到達時は早期終了する。
    反応候補がなくても両MDを続け、次cycleで近接を再探索する。
    開始点数／割合は引数を省略するとactivation.yamlを使用する。
    max_reactions_per_cycleを指定するとYAMLのルール上限もその値へ上書きする。
    reaction_cutoff: Å。YAMLの距離を上書きし、既定5 Åから始める。
    反応0件がcutoff_zero_cycles回連続すると次cycleからcutoff_growth_factor倍。
    max_reaction_cutoffで頭打ち。固定距離ならfixed_reaction_cutoff=True。
    遠い対を無理に結合させると歪んだ構造になりやすいので、長いrunでは上限を5〜6 Åに下げる。
    relaxation_stages: 各cycleの反応後MDの段リスト（cycle_stages参照）。Noneならnvt_*/npt_*の2段。
    stop_after_idle_cycles: 探索距離が上限（固定距離ならその値）のまま反応0件のcycleが
      この回数続いたら終了する（既定200、0で無効）。活性点が網目に埋もれて動けなくなった後の空回りを防ぐ。
    reaction_minimization_radius: 反応後の最小化で動かす範囲 Å（既定12）。反応原子からこの距離内と
      更新した結合近傍だけを最小化し、他は固定する。0なら系全体を最小化する（系全体の熱構造も消える）。
    pack_relax / pack_postrelax: packing中と後のnativeの局所緩和（PAVES既定と同じ50 sweep、400）。
      0にすると生の成長座標が残り、結合が環を貫通（ring spearing）して最小化でもほどけない。
      packing後に環を貫通する結合が1本でもあれば停止する。
    final_relaxation_steps / final_md_steps: 反応終了後、最終構造を保存する前に
      NVT 0.1 fs × final_relaxation_steps（既定1000）→ NPT 1 atm・1 fs × final_md_steps（既定1,000,000 = 1 ns）。
      反応を終えた直後の構造はreacted.structure.json、緩和後をfinal.*として保存する。0で各段を省略。
    template_workers: 未知の反応テンプレートを並列計算するプロセス数（既定4、1で逐次）。
      結果は逐次計算と同一。1プロセスあたり約1 GBのメモリを使う。
    post_reaction_md_steps: 反応直後の最小化に続けて、反応の確定前にnvt_timestep刻みで
      流すNVT step数（既定0）。commit前に伸びた新結合をならす。ここで非有限エネルギーに
      なると反応はcommitされず、計算はReactionTransactionErrorで止まる。
    reaction_assignment: local_templatesは局所電荷転写／パラメータ更新、wholeは全体再計算。
    local_radius: 更新対象の結合本数。local_context_radius: 計算用の追加外側buffer。
    template_cache: 計算済みフラグメントbank。省略時results/<project>/reaction_templates。
    既知の環境は再計算せず、未知の環境は上限原子数以内で一度計算して保存。
    require_reaction_templates=Trueなら、成長cycle中の未知templateを計算せずエラー。
    target_conversionは0〜1（既定0.995）、Noneならcyclesまで。
    save_trajectory: TrueならOVITO軌跡とcycle別構造も保存する（既定False）。
    remove_sol_fraction: Trueなら現像を模して、周期境界をまたいで無限につながった網目だけを
      最終構造に残し、未反応モノマーとオリゴマー（ゾル分）を削除する（既定False）。
      削除前はfinal_with_sol.*に残す。抜けた体積はこのスクリプトでは緩和しない。
    最終構造と反応履歴は常に保存。集計JSONはcycleごとに更新する。
    戻り値はReactiveMDResult。
    """
    if not isinstance(project, str) or not re.fullmatch(r'[\w][\w.-]*', project):
        raise ValueError('project must be a name using letters, digits, underscore, dot or hyphen; no paths')
    if type(force) is not bool:
        raise ValueError('force must be a boolean')
    if type(save_trajectory) is not bool:
        raise ValueError('save_trajectory must be a boolean')
    if type(fixed_reaction_cutoff) is not bool:
        raise ValueError('fixed_reaction_cutoff must be a boolean')
    cutoff_policy = AdaptiveReactionCutoff(initial_cutoff=reaction_cutoff,
        maximum_cutoff=max_reaction_cutoff, growth_factor=cutoff_growth_factor,
        zero_reaction_cycles=cutoff_zero_cycles)
    definitions = Path.cwd()
    results_directory = definitions / 'results'
    output = results_directory / project
    if results_directory.is_symlink() or output.is_symlink():
        raise ValueError('results and project directories must not be symbolic links')
    if output.exists() and not force:
        raise FileExistsError(f'Project already exists: {output}. Use --force to overwrite it.')
    if reaction_assignment not in ('local_templates', 'whole'):
        raise ValueError('reaction_assignment must be local_templates or whole')
    bank = output / 'reaction_templates' if template_cache is None else Path(template_cache).absolute()
    field = (LocalSageForceField(p.SageForceField(version='2.3.0'), radius=local_radius,
             context_radius=local_context_radius, max_fragment_atoms=max_fragment_atoms, cache_dir=bank,
             require_templates=require_reaction_templates, workers=template_workers)
             if reaction_assignment == 'local_templates' else p.SageForceField(version='2.3.0'))

    from rdkit import Chem
    from polypaves.config import SimulationConfig
    from polypaves.api import _number

    if type(total_atoms) is not int or total_atoms < 1:
        raise ValueError('total_atoms must be a positive integer')
    if activation_count is not None and activation_fraction is not None:
        raise ValueError('Use activation_count or activation_fraction, not both')
    if max_reactions_per_cycle is not None and (type(max_reactions_per_cycle) is not int or max_reactions_per_cycle < 1):
        raise ValueError('max_reactions_per_cycle must be a positive integer')
    if atom_count_mode not in ('above', 'nearest', 'floor', 'exact'):
        raise ValueError('atom_count_mode must be above, nearest, floor or exact')
    SimulationConfig(temperature=temperature, timestep=nvt_timestep,
                     nvt_steps=nvt_steps, npt_steps=npt_steps)
    _number(pressure_atm, 'pressure_atm')
    _number(density, 'density')
    _number(reaction_minimization_radius, 'reaction_minimization_radius')
    if reaction_minimization_radius < 0:
        raise ValueError('reaction_minimization_radius must be >= 0')
    for name, value, minimum in (('cycles', cycles, 0),
                                 ('minimization_iterations', minimization_iterations, 0),
                                 ('compression_steps', compression_steps, 0),
                                 ('post_reaction_md_steps', post_reaction_md_steps, 0),
                                 ('stop_after_idle_cycles', stop_after_idle_cycles, 0),
                                 ('template_workers', template_workers, 1),
                                 ('final_relaxation_steps', final_relaxation_steps, 0),
                                 ('final_md_steps', final_md_steps, 0),
                                 ('pack_relax', pack_relax, 0), ('pack_postrelax', pack_postrelax, 0)):
        if type(value) is not int or value < minimum:
            raise ValueError(f'{name} must be an integer >= {minimum}')
    if target_conversion is not None and not 0 <= target_conversion <= 1:
        raise ValueError('target_conversion must be between 0 and 1')
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError('smiles must describe one connected epoxy molecule')
    epoxides = molecule.GetSubstructMatches(Chem.MolFromSmarts('[C]1[O;+0;X2][C]1'))
    if not epoxides:
        raise ValueError('smiles must contain a neutral C-O-C epoxide ring')
    # Kekulé入力（C1=CC=CC=C1など）はそのまま単結合/二重結合のgraphになる。反応後の
    # 再検証で原子だけが芳香族になり、局所フラグメントが環の単結合で切れて
    # "non-ring atom N marked aromatic"で停止するため、芳香族表記の正規SMILESでpackする。
    input_smiles, smiles = smiles, Chem.MolToSmiles(molecule)
    atoms_per_molecule = Chem.AddHs(molecule).GetNumAtoms()
    activation_rule = load_activation_definitions(definitions / 'activation.yaml')[0]
    rules = load_reaction_definitions(definitions / 'reactions.yaml')
    selector = rules.selection
    cli_selection = dict(count=activation_count, fraction=activation_fraction)
    if activation_count is None and activation_fraction is None:
        if activation_rule.selection in ('fraction', 'random') and activation_rule.fraction is not None:
            activation_fraction = activation_rule.fraction
        elif activation_rule.selection in ('count', 'random'):
            activation_count = activation_rule.count
        else:
            raise ValueError('This example supports count or fraction activation selection')
    sites_per_molecule = len(Chem.AddHs(molecule).GetSubstructMatches(Chem.MolFromSmarts(activation_rule.target_smarts)))
    plan = plan_size(total_atoms, atoms_per_molecule, sites_per_molecule,
                     added_per_site=len(activation_rule.atoms_to_add), count=activation_count,
                     fraction=activation_fraction, mode=atom_count_mode)
    if max_reactions_per_cycle is not None:
        rules = [replace(rule, max_reactions_per_cycle=max_reactions_per_cycle) for rule in rules]
    rules = [replace(rule, distance_cutoff=reaction_cutoff) for rule in rules]
    pressure_bar = pressure_atm * 1.01325
    if relaxation_stages is not None:
        stages = cycle_stages(relaxation_stages, temperature=temperature, pressure_atm=pressure_atm)
    else:
        stages = (MDStage('NVT', steps=nvt_steps, timestep=nvt_timestep, temperature=temperature),
                  MDStage('NPT', steps=npt_steps, timestep=npt_timestep, temperature=temperature,
                          pressure=pressure_bar))
    if compression_steps:
        _number(compression_pressure_mpa, 'compression_pressure_mpa')
        if compression_pressure_mpa <= 0:
            raise ValueError('compression_pressure_mpa must be positive')
    # 1 MPa = 10 bar。加圧後は同じ座標・速度・boxのまま1 atmへ戻す。
    compression = (MDStage('NPT', steps=compression_steps, timestep=equilibration_timestep,
                           temperature=temperature, pressure=compression_pressure_mpa * 10.)
                   if compression_steps else None)
    equilibration = MDStage('NPT', steps=equilibration_steps, timestep=equilibration_timestep,
                            temperature=temperature, pressure=pressure_bar)

    # 入力と両YAMLを検証してから、指定されたprojectだけを作り直す。
    if results_directory.is_symlink() or output.is_symlink():
        raise ValueError('results and project directories must not be symbolic links')
    if output.exists():
        if not force:
            raise FileExistsError(f'Project already exists: {output}. Use --force to overwrite it.')
        if not output.is_dir():
            raise NotADirectoryError(f'Project path is not a directory: {output}')
        shutil.rmtree(output)
    output.mkdir(parents=True)
    shutil.copy2(definitions / 'reactions.yaml', output / 'reactions.yaml')
    # CLIの開始点数／割合を保存先のactivation.yamlにも反映し、実際の設定と一致させる。
    (output / 'activation.yaml').write_text(activation_yaml_with_selection(
        (definitions / 'activation.yaml').read_text(), **cli_selection))
    precision = 'mixed' if platform == 'CUDA' else None

    # p.Solventはpackへの分子指定。開始点で追加するHを原子予算から引く。
    print(f"Project {project}: {plan['monomers']} molecules, {plan['actual_total_atoms']} atoms, "
          f"{plan['activation_count']}/{plan['initial_epoxide_sites']} epoxide starts", flush=True)
    structure = p.pack(
        [p.Solvent(smiles, name='epoxy')], counts=[plan['monomers']],
        density=density, temperature=temperature,
        forcefield='sage', forcefield_version='2.3.0',
        platform=platform, precision=precision, seed=seed,
        trials=20, relax=pack_relax, postrelax=pack_postrelax,
    )
    packed_atoms = structure.n_atoms
    structure.write_json(output / 'before_activation.structure.json')
    piercing = structure.report.get('ring_piercing_bonds')
    if piercing is None:
        raise RuntimeError('This polypaves does not report ring piercing after packing')
    if piercing:
        raise RuntimeError(f'{piercing} bonds pierce rings after packing; increase --pack-postrelax/--pack-relax '
                           'or lower --density (see before_activation.structure.json)')
    if packed_atoms != plan['packed_atoms']:
        raise RuntimeError('Packed atom count differs from the size plan')
    forcefield = field.base if isinstance(field, LocalSageForceField) else field
    initializer = ActivationInitializer(structure, forcefield=forcefield, seed=seed)
    bridge = OpenMMBridge(platform=platform, precision=precision,
                          temperature=temperature, timestep=nvt_timestep, seed=seed)

    def validate_initialization(graph, activation):
        return bridge.build(activation.parameterized, graph)

    selection = {'fraction': activation_fraction} if activation_fraction is not None else {'count': activation_count}
    activated = initializer.apply(activation_rule, **selection, validator=validate_initialization)
    if structure.n_atoms != plan['actual_total_atoms'] or activated.number_of_activations != plan['activation_count']:
        raise RuntimeError('Activation atom count differs from the packing budget')
    initializer.history.write_jsonl(output / 'activation_history.jsonl')
    runner = ReactiveMDRunner(
        structure=structure, forcefield=field, reactions=rules,
        selector=selector, temperature=temperature, timestep=nvt_timestep,
        platform=platform, precision=precision, seed=seed,
        adaptive_cutoff=None if fixed_reaction_cutoff else cutoff_policy,
        dirty_radius=local_radius if reaction_assignment == 'local_templates' else 4,
    )
    runner.prepare(minimize=True, minimization_max_iterations=minimization_iterations)
    if isinstance(field, LocalSageForceField):
        # 初期構造の反応候補から、反応後の小フラグメントをMD前に準備する。
        # 現在のgraphは変更せず、未知の環境も全体計算にはfallbackしない。
        selected = copy.deepcopy(runner.engine.selector).select(runner.engine.find_candidates(),
                                        max_reactions=max_reactions_per_cycle)
        field.precompute(runner.engine.graph, selected)
    runner.structure.write_json(output / 'initial.structure.json')
    # 開始点の付加後、反応前の分子番号。削除したゾル分のモノマー／オリゴマーの区別に使う。
    origin = {a['stable_atom_id']: i for i, g in enumerate(runner.structure.graphs) for a in g['atoms']}
    settings = dict(project=project, working_directory=str(definitions), output_directory=str(output),
                    smiles=smiles, input_smiles=input_smiles, target_total_atoms=total_atoms,
                    atom_count_mode=atom_count_mode, atoms_per_molecule=atoms_per_molecule,
                    packed_atoms=packed_atoms, actual_total_atoms=structure.n_atoms,
                    monomers=plan['monomers'], activation_count=plan['activation_count'],
                    activation_fraction=activation_fraction, initial_epoxide_sites=plan['initial_epoxide_sites'],
                    actual_activation_fraction=plan['activation_count'] / plan['initial_epoxide_sites'],
                    density_g_cm3=density, temperature_K=temperature, pressure_atm=pressure_atm,
                    compression=asdict(compression) if compression else None,
                    compression_pressure_MPa=compression_pressure_mpa if compression else None,
                    equilibration=asdict(equilibration),
                    cycle_md_stages=[asdict(stage) for stage in stages], cycles=cycles,
                    relaxation_stages=relaxation_stages, post_reaction_md_steps=post_reaction_md_steps,
                    stop_after_idle_cycles=stop_after_idle_cycles, template_workers=template_workers,
                    reaction_minimization_radius_A=reaction_minimization_radius or None,
                    final_relaxation_steps=final_relaxation_steps, final_md_steps=final_md_steps,
                    pack_relax=pack_relax, pack_postrelax=pack_postrelax,
                    ring_piercing_bonds_after_packing=piercing,
                    save_trajectory=save_trajectory, remove_sol_fraction=remove_sol_fraction,
                    output_interval_cycles=1 if save_trajectory else None,
                    minimization_iterations=minimization_iterations,
                    target_conversion=target_conversion, max_reactions_per_cycle=max_reactions_per_cycle,
                    adaptive_cutoff=None if fixed_reaction_cutoff else asdict(cutoff_policy),
                    reaction_assignment=field.statistics if isinstance(field, LocalSageForceField) else {'mode': 'whole'},
                    fragment_charge_transfer_approximation=isinstance(field, LocalSageForceField),
                    initial_distance_cutoffs_A={rule.name: rule.distance_cutoff for rule in runner.engine.reactions},
                    effective_reaction_limits={rule.name: rule.max_reactions_per_cycle for rule in rules}, seed=seed)
    (output / 'settings.json').write_text(json.dumps(settings, indent=2) + '\n')

    with (OVITOTrajectoryWriter(output, runner.engine.graph)
          if save_trajectory else nullcontext()) as trajectory:
        def reactive_atom_ids():
            return {sid for roles in runner.engine.tracker.current.values()
                    for sites in roles for site in sites for sid in site}

        def write_frame(active, **metadata):
            if trajectory is not None:
                trajectory.write(active.engine.graph,
                    OpenMMStateSnapshot.from_simulation(active.simulation),
                    reactive_atom_ids=reactive_atom_ids(), **metadata)
                trajectory.flush()

        def save_summary(active):
            metrics = {**active.metrics, 'frames': trajectory.num_frames if trajectory is not None else 0,
                       'initiated_sites': plan['activation_count'],
                       'initial_epoxide_sites': plan['initial_epoxide_sites'],
                       'epoxide_consumption_including_activation':
                           (plan['activation_count'] + len(active.history)) / plan['initial_epoxide_sites'],
                       'target_total_atoms': total_atoms, 'actual_total_atoms': structure.n_atoms,
                       'atom_count_difference': structure.n_atoms - total_atoms,
                       'platform': active.simulation.context.getPlatform().getName()}
            if isinstance(field, LocalSageForceField):
                metrics['reaction_templates'] = field.statistics
            (output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
            (output / 'cycles.json').write_text(json.dumps(active.cycle_history, indent=2) + '\n')
            return metrics

        def save_structure(active, prefix):
            active.structure.write_json(Path(str(prefix) + '.structure.json'))
            write_structure_files(active.structure, prefix, charges=[
                atom['partial_charge'] for graph in active.structure.graphs for atom in graph['atoms']])

        write_frame(runner, phase='initial')
        # 初期緩和中は成長反応を実行しない。座標・速度・boxを成長へ引き継ぐ。
        runner.equilibrate([compression, equilibration] if compression else [equilibration])
        runner.structure.write_json(output / 'before_reaction.structure.json')
        write_frame(runner, phase='equilibrated')
        save_summary(runner)

        def save_cycle(active, cycle):
            write_frame(active,
                cycle=cycle['cycle'], phase='cycle end', reactions=len(active.history),
                conversion=cycle['conversion'])
            if save_trajectory:
                # 有効時のみ、cycleごとの構造も可視化用の履歴として保存する。
                directory = output / 'cycles'
                directory.mkdir(exist_ok=True)
                save_structure(active, directory / f"cycle_{cycle['cycle'] + 1:05d}")
            save_summary(active)

        # md_steps_per_cycle=0で各cycleを反応探索から始める。
        # 軌跡と構造の出力はon_cycleでまとめて行う。
        result = runner.run(md_steps_per_cycle=0, max_cycles=cycles,
            max_reactions_per_cycle=max_reactions_per_cycle, minimize_after_reaction=True,
            minimization_max_iterations=minimization_iterations,
            post_reaction_md_steps=post_reaction_md_steps, cycle_md_stages=stages,
            stop_after_idle_cycles=stop_after_idle_cycles or None,
            minimization_radius=reaction_minimization_radius or None, target_conversion=target_conversion,
            history_path=output / 'reaction_history.jsonl', on_cycle=save_cycle)
        save_summary(runner)
        # 反応直後の構造を残し、0.1 fsの緩和と1 atmのMDを経た構造を最終構造とする。
        runner.structure.write_json(output / 'reacted.structure.json')
        final_stages = []
        if final_relaxation_steps:
            final_stages.append(MDStage('NVT', steps=final_relaxation_steps, timestep=.1, temperature=temperature))
        if final_md_steps:
            final_stages.append(MDStage('NPT', steps=final_md_steps, timestep=1., temperature=temperature,
                                        pressure=pressure_bar))
        if final_stages:
            runner.equilibrate(final_stages)
            write_frame(runner, phase='final relaxed')
        metrics = save_summary(runner)
        final = runner.structure
        if remove_sol_fraction:
            # 削除前の構造も残す。抜けた体積は後続のNPT（network_propertiesの熱履歴など）で緩和する。
            save_structure(runner, output / 'final_with_sol')
            final, sol = without_sol(final, origin)
            metrics.update(sol)
            (output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
            if final is None:
                raise RuntimeError('No molecule percolates through the cell (below the gel point): '
                                   'the whole film would dissolve. See final_with_sol.structure.json')
            print(f"Removed the sol fraction: {sol['sol_monomers_removed']} monomers and "
                  f"{sol['sol_oligomers_removed']} oligomers ({sol['sol_atoms_removed']} atoms, "
                  f"{100 * sol['sol_mass_fraction']:.1f} wt%); {final.n_atoms} atoms remain", flush=True)
        save_structure(SimpleNamespace(structure=final), output / 'final')
        # 環を貫通した結合は伸びたまま残り、後続の固定トポロジーMDの結合長検査で停止する。
        from polypaves.ring_piercing import ring_piercing_bonds
        threaded = ring_piercing_bonds(final.graphs, final._positions, final.box_vectors)
        metrics['ring_piercing_bonds_final'] = len(threaded)
        (output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
        if threaded:
            raise RuntimeError(f'{len(threaded)} bonds pierce rings in the final structure '
                               f'(atom indices {threaded[:5]}); it cannot be used for later MD')
    if isinstance(field, LocalSageForceField):
        field.close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument('--smiles', help='エポキシ1分子のSMILES')
    inputs.add_argument('--smiles-file', type=Path, help='SMILESを1つだけ書いたテキストファイル')
    parser.add_argument('--total-atoms', type=int, default=DEFAULT_TOTAL_ATOMS,
                        help='開始点の追加Hを含む目標総原子数（既定20000）')
    parser.add_argument('--atom-count-mode', choices=('above', 'nearest', 'floor', 'exact'),
                        default=DEFAULT_ATOM_COUNT_MODE,
                        help='原子数の丸め方。aboveは目標を超える最小分子数（既定above）')
    parser.add_argument('--density', type=float, default=.5,
                        help='packing時の密度 g/cm^3（既定0.5。加圧NPTで縮める）')
    starts = parser.add_mutually_exclusive_group()
    starts.add_argument('--activation-count', type=int)
    starts.add_argument('--activation-fraction', type=float, help='開始点の割合。省略時はactivation.yamlを使用')
    parser.add_argument('--project', required=True, help='保存先: カレントのresults/<project>')
    parser.add_argument('--force', action='store_true', help='既存projectの内容を削除して再実行する')
    parser.add_argument('--save-trajectory', action='store_true',
                        help='OVITO軌跡とcycle別構造を保存する（既定OFF）')
    parser.add_argument('--platform', choices=('CUDA', 'CPU'), default='CUDA')
    parser.add_argument('--cycles', type=int, default=DEFAULT_CYCLES,
                        help='最大反応サイクル数（既定5000）')
    parser.add_argument('--temperature', type=float, default=300.)
    parser.add_argument('--compression-steps', type=int, default=50_000,
                        help='packing後の加圧NPTのstep数（既定50000）。既定1 fs刻みで50 ps、0で無効')
    parser.add_argument('--compression-pressure-mpa', type=float, default=200.,
                        help='加圧NPTの圧力 MPa（既定200）')
    parser.add_argument('--equilibration-steps', type=int, default=1_000_000,
                        help='初期NPTのstep数（既定1000000）。既定1 fs刻みで1 ns')
    parser.add_argument('--equilibration-timestep', type=float, default=1.)
    parser.add_argument('--nvt-timestep', type=float, default=.1)
    parser.add_argument('--nvt-steps', type=int, default=10000)
    parser.add_argument('--npt-timestep', type=float, default=1.)
    parser.add_argument('--npt-steps', type=int, default=1000)
    parser.add_argument('--pressure-atm', type=float, default=1.)
    parser.add_argument('--minimization-iterations', type=int, default=2000)
    parser.add_argument('--target-conversion', type=float, default=DEFAULT_TARGET_CONVERSION,
                        help='目標反応率0〜1（既定0.995 = 99.5%%）')
    parser.add_argument('--max-reactions-per-cycle', type=int, help='1 cycleの最大反応数。指定時はYAMLのルール上限も上書き')
    parser.add_argument('--reaction-cutoff', type=float,
                        help='開始時の反応探索距離 Å（YAMLを上書き、既定5。上限より大きければ上限）')
    parser.add_argument('--max-reaction-cutoff', type=float, default=7., help='探索距離の上限 Å（既定7）')
    parser.add_argument('--cutoff-growth-factor', type=float, default=1.2, help='0件が続いた際の距離倍率（既定1.2）')
    parser.add_argument('--cutoff-zero-cycles', type=int, default=2, help='距離拡大までの連続0件cycle数（既定2）')
    parser.add_argument('--fixed-reaction-cutoff', action='store_true', help='距離をreaction-cutoffに固定し、自動拡大を無効化')
    parser.add_argument('--relaxation-stages', type=load_relaxation_stages,
                        help='各cycleの反応後MD。JSON文字列かJSON/YAMLファイル。指定時は--nvt-*/--npt-*を使わない')
    parser.add_argument('--post-reaction-md-steps', type=int, default=0,
                        help='反応確定前、最小化の後に流すNVT step数（--nvt-timestep刻み、既定0）')
    parser.add_argument('--stop-after-idle-cycles', type=int, default=200,
                        help='探索距離が上限のまま反応0件がこのcycle数続いたら終了（既定200、0で無効）')
    parser.add_argument('--reaction-minimization-radius', type=float, default=12.,
                        help='反応後の最小化で動かす範囲 Å（既定12、0で系全体）')
    parser.add_argument('--pack-relax', type=int, default=50, help='packing中の局所緩和sweep数（既定50）')
    parser.add_argument('--pack-postrelax', type=int, default=400, help='packing後の局所緩和（既定400）')
    parser.add_argument('--final-relaxation-steps', type=int, default=1000,
                        help='反応終了後のNVT 0.1 fsのstep数（既定1000、0で省略）')
    parser.add_argument('--final-md-steps', type=int, default=1_000_000,
                        help='その後のNPT 1 atm・1 fsのstep数（既定1000000 = 1 ns、0で省略）')
    parser.add_argument('--template-workers', type=int, default=4,
                        help='未知の反応テンプレートを並列計算するプロセス数（既定4、1で逐次）')
    parser.add_argument('--reaction-assignment', choices=('local_templates', 'whole'), default='local_templates',
                        help='反応後の割り当て方式。既定は局所フラグメント電荷の転写')
    parser.add_argument('--local-radius', type=int, default=5, help='更新対象の結合本数（既定5）')
    parser.add_argument('--local-context-radius', type=int, default=8, help='更新領域の外側に含める結合本数（既定8）')
    parser.add_argument('--max-fragment-atoms', type=int, default=1200, help='1フラグメントの計算原子数上限（既定1200）')
    parser.add_argument('--template-cache', type=Path, help='計算済み電荷／パラメータbank。省略時project内reaction_templates')
    parser.add_argument('--require-reaction-templates', action='store_true', help='成長cycle中の未知templateを計算せずエラーにする')
    parser.add_argument('--seed', type=int, default=31415)
    parser.add_argument('--remove-sol-fraction', action='store_true',
                        help='現像を模して、周期境界をまたいでつながった網目だけを残し、未反応モノマーとオリゴマーを削除する')
    parser.add_argument('--log-interval', type=float, default=10.)
    args = vars(parser.parse_args())
    smiles_file = args.pop('smiles_file')
    args['smiles'] = smiles_file.read_text().strip() if smiles_file else args['smiles'] or DEFAULT_SMILES
    if args['reaction_cutoff'] is None:
        args['reaction_cutoff'] = min(5., args['max_reaction_cutoff'])
    p.configure_logging(progress_interval=args.pop('log_interval'))
    run(**args)


if __name__ == '__main__':
    main()
