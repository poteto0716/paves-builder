"""Equilibrium statistics of the unreacted reactive rings (epoxide, oxetane, ...) before any reaction.

Each unreacted ring is one group: its ring oxygen, the attackable carbons
(ring carbons bonded to that oxygen carrying two H, the CH2 the relay
reaction opens) and its ring heavy atoms. Rings come either from the
reaction tracker (site_groups: the remaining sites of a rule's conversion
reactant) or from the activation rule's target_smarts matched on plain
chemical graphs (smarts_groups: an unreacted monomer liquid, no reaction
engine). The definition does not depend on the ring size: a three-membered
ring is reported as epoxide, a four-membered one as oxetane.

Two groups are in contact when the oxygen of either is within the cutoff of
an attackable carbon of the other: the geometry of one relay step (the OH
left on an opened ring's oxygen attacks the next ring's CH2; for oxetane the
next ring's oxygen attacks the alpha-CH2 of the oxonium, the same pair).
Intramolecular pairs count as well. Distances use the minimum image.

Per frame: number density (groups / cell volume), coordination number (mean
contacts per group), largest cluster fraction of the contact graph (largest
connected group count / all groups). MSD uses the ring heavy-atom centroid,
unwrapped frame to frame, with the drift of the system centre of mass removed
and averaged over all time origins.
"""
import csv
import json
from pathlib import Path

import numpy as np

A2_PS_TO_M2_S = 1e-8
PER_NM3_TO_MMOL_CM3 = 1 / .60221408  # 1 nm^-3 = 1/0.6022 mmol/cm^3
RING_NAMES = {3: 'epoxide', 4: 'oxetane'}


def ring_name(size):
    return RING_NAMES.get(size, f'{size}-membered ring')


def _group(heavy, atomic_number, neighbors, hydrogens):
    """Group of one ring from its heavy atoms; None unless it has one O and an attackable CH2."""
    oxygens = [a for a in heavy if atomic_number(a) == 8]
    if len(oxygens) != 1:
        return None
    oxygen = oxygens[0]
    carbons = [a for a in heavy if atomic_number(a) == 6 and oxygen in neighbors(a) and hydrogens(a) == 2]
    return dict(oxygen=oxygen, carbons=carbons, heavy=list(heavy)) if carbons else None


def site_groups(graph, sites, order):
    """Index sets (into the simulation atom order) of each unreacted ring of the reaction tracker."""
    index = {sid: i for i, sid in enumerate(order)}
    z = lambda sid: graph.atoms[sid]['atomic_number']
    groups = []
    for site in sorted(sites):
        heavy = [sid for sid in site if z(sid) != 1]
        group = _group(heavy, z, lambda sid: graph.neighbors[sid],
                       lambda sid: sum(z(j) == 1 for j in graph.neighbors[sid]))
        if group:
            groups.append(dict(oxygen=index[group['oxygen']], carbons=[index[c] for c in group['carbons']],
                               heavy=[index[sid] for sid in heavy]))
    return groups


def smarts_groups(graphs, target_smarts):
    """Rings matching target_smarts in plain graphs, indexed in the concatenated atom order.

    Each ring is counted once whatever the number of SMARTS mappings onto it
    (an oxetane with two ring CH2 matches twice).
    """
    from rdkit import Chem
    from polypaves.molecular import _rdkit_molecule
    query = Chem.MolFromSmarts(target_smarts)
    if query is None:
        raise ValueError(f'invalid target_smarts: {target_smarts}')
    groups, offset, cache = [], 0, {}
    for graph in graphs:
        atoms, bonds = graph['atoms'], graph['bonds']
        # Identical monomers share one match; keyed by the chemistry, not the atom IDs.
        key = (tuple((a['atomic_number'], a['formal_charge'], a['aromatic']) for a in atoms),
               tuple(tuple(b) for b in bonds))
        if key not in cache:
            neighbors = [set() for _ in atoms]
            for i, j, _ in bonds:
                neighbors[i].add(j)
                neighbors[j].add(i)
            z = lambda i: atoms[i]['atomic_number']
            molecule = _rdkit_molecule(graph, perceive=False)
            found, seen = [], set()
            for match in molecule.GetSubstructMatches(query, uniquify=True, maxMatches=100000):
                heavy = sorted(i for i in match if z(i) != 1)
                if frozenset(heavy) in seen:
                    continue
                seen.add(frozenset(heavy))
                group = _group(heavy, z, lambda i: neighbors[i], lambda i: sum(z(j) == 1 for j in neighbors[i]))
                if group:
                    found.append(group)
            cache[key] = found
        for group in cache[key]:
            groups.append(dict(oxygen=group['oxygen'] + offset, carbons=[c + offset for c in group['carbons']],
                               heavy=[a + offset for a in group['heavy']]))
        offset += len(atoms)
    return groups


def _minimum_image(delta, box, inverse):
    fractional = delta @ inverse
    return (fractional - np.rint(fractional)) @ box


def contacts(oxygen, carbons, box, cutoff, chunk=256):
    """Boolean contact matrix: O_i within cutoff of an attackable C_j, or O_j of C_i."""
    inverse = np.linalg.inv(box)
    n = len(oxygen)
    near = np.zeros((n, n), dtype=bool)
    for slot in range(carbons.shape[1]):
        column = carbons[:, slot]
        valid = np.isfinite(column[:, 0])
        for start in range(0, n, chunk):
            delta = column[None, :, :] - oxygen[start:start + chunk, None, :]
            distance = np.linalg.norm(_minimum_image(delta, box, inverse), axis=-1)
            near[start:start + chunk] |= (distance <= cutoff) & valid[None, :]
    near |= near.T
    np.fill_diagonal(near, False)
    return near


def cluster_statistics(near):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    count, labels = connected_components(csr_matrix(near), directed=False)
    sizes = np.bincount(labels)
    return dict(coordination_number=float(near.sum(axis=1).mean()),
                largest_cluster_fraction=float(sizes.max() / len(labels)), clusters=int(count))


def mean_squared_displacement(centroids, times_ps):
    """MSD(lag) averaged over all time origins; centroids is frames x groups x 3 (angstrom)."""
    frames = len(centroids)
    lags = np.arange(1, frames)
    msd = np.array([np.mean(np.sum((centroids[k:] - centroids[:-k]) ** 2, axis=-1)) for k in lags])
    lag_ps = np.asarray(times_ps)[lags] - times_ps[0]
    return lag_ps, msd


def diffusion_fit(lag_ps, msd):
    """Einstein slope over lags 10-50 % of the run: D (m2/s) and the log-log exponent alpha."""
    if not len(lag_ps):
        return dict(diffusion_coefficient_m2_s=None, msd_exponent_alpha=None, fit_window_ps=None)
    window = (lag_ps >= .1 * lag_ps[-1]) & (lag_ps <= .5 * lag_ps[-1]) & (msd > 0)
    if window.sum() < 3:
        return dict(diffusion_coefficient_m2_s=None, msd_exponent_alpha=None, fit_window_ps=None)
    slope = np.polyfit(lag_ps[window], msd[window], 1)[0]
    alpha = np.polyfit(np.log(lag_ps[window]), np.log(msd[window]), 1)[0]
    return dict(diffusion_coefficient_m2_s=float(slope / 6 * A2_PS_TO_M2_S), msd_exponent_alpha=float(alpha),
                fit_window_ps=[float(lag_ps[window][0]), float(lag_ps[window][-1])])


class RingSampler:
    """Collects frames (wrapped or not, angstrom) and unwraps them frame to frame."""

    def __init__(self, groups, masses, cutoffs):
        self.groups, self.cutoffs = groups, sorted(set(float(c) for c in cutoffs))
        self.masses = np.asarray(masses, dtype=float)
        width = max(len(g['carbons']) for g in groups)
        self.carbon_index = np.array([g['carbons'] + [-1] * (width - len(g['carbons'])) for g in groups])
        self.oxygen_index = np.array([g['oxygen'] for g in groups])
        self.rows, self.centroids, self.times = [], [], []
        self.previous = self.unwrapped = None

    def add(self, positions, box, time_ps):
        positions, box = np.asarray(positions, dtype=float), np.asarray(box, dtype=float)
        inverse = np.linalg.inv(box)
        if self.previous is None:
            self.unwrapped = positions.copy()
            # Whole rings at the first frame: each atom at its minimum image from the ring's first atom.
            self.start = positions.copy()
            self.centroid0 = np.array([positions[g['heavy'][0]] + _minimum_image(
                positions[g['heavy']] - positions[g['heavy'][0]], box, inverse).mean(axis=0) for g in self.groups])
        else:
            self.unwrapped += _minimum_image(positions - self.previous, box, inverse)
        self.previous = positions
        moved = self.unwrapped - self.start
        drift = self.masses @ moved / self.masses.sum()
        self.centroids.append(self.centroid0 + np.array([moved[g['heavy']].mean(axis=0) for g in self.groups]) - drift)
        self.times.append(float(time_ps))
        oxygen = positions[self.oxygen_index]
        carbons = np.where(self.carbon_index[..., None] >= 0, positions[self.carbon_index], np.nan)
        volume = abs(float(np.linalg.det(box))) / 1000.
        row = dict(time_ps=float(time_ps), volume_nm3=volume, rings=len(self.groups),
                   ring_number_density_nm3=len(self.groups) / volume)
        for cutoff in self.cutoffs:
            stats = cluster_statistics(contacts(oxygen, carbons, box, cutoff))
            row.update({f'{key}_{cutoff:g}A': value for key, value in stats.items()})
        self.rows.append(row)

    def summary(self, output, *, temperature_K, pressure_atm, sampling_ps):
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        lag_ps, msd = mean_squared_displacement(np.array(self.centroids), self.times)

        def statistics(key):
            values = np.array([row[key] for row in self.rows])
            return dict(mean=float(values.mean()), std=float(values.std()))

        by_cutoff = {f'{cutoff:g}': {key: statistics(f'{key}_{cutoff:g}A')
                                     for key in ('coordination_number', 'largest_cluster_fraction', 'clusters')}
                     for cutoff in self.cutoffs}
        size = len(self.groups[0]['heavy'])
        density = statistics('ring_number_density_nm3')
        result = dict(temperature_K=temperature_K, pressure_atm=pressure_atm, sampling_ps=sampling_ps,
                      frames=len(self.rows), rings=len(self.groups), ring_size=size, ring=ring_name(size),
                      volume_nm3=statistics('volume_nm3'), ring_number_density_nm3=density,
                      ring_number_density_mmol_cm3={k: v * PER_NM3_TO_MMOL_CM3 for k, v in density.items()},
                      contact_cutoffs_A=self.cutoffs, by_cutoff_A=by_cutoff,
                      msd_final_A2=float(msd[-1]) if len(msd) else None,
                      msd_final_lag_ps=float(lag_ps[-1]) if len(lag_ps) else None,
                      **diffusion_fit(lag_ps, msd),
                      definitions='group = unreacted reactive ring (activation target); contact = ring O of one '
                                  'within the cutoff of an attackable CH2 carbon of the other (minimum image, '
                                  'intramolecular included); coordination = mean contacts per group; largest '
                                  'cluster fraction = largest connected group count / groups; MSD of the ring '
                                  'heavy-atom centroid, centre-of-mass drift removed, all time origins; D from '
                                  'MSD/6t over lags 10-50 % of the run (check msd_exponent_alpha ~ 1).')
        (output / 'reactive_rings.json').write_text(json.dumps(result, indent=2) + '\n')
        with open(output / 'reactive_ring_frames.csv', 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(self.rows[0]))
            writer.writeheader()
            writer.writerows(self.rows)
        with open(output / 'reactive_ring_msd.csv', 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['lag_ps', 'msd_A2'])
            writer.writerows(zip(lag_ps.tolist(), msd.tolist()))
        return result
