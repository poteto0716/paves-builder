"""Elastically effective strands, junctions, dangling material and polar groups of a network.

The analysis works on the heavy-atom graph (hydrogen mass is carried by its
heavy atom) with each small-ring system (fused 3- to 8-membered rings)
contracted to one node, so an aromatic core in a chain is a pass-through,
not a set of branch points.

1. Gel: molecules that wrap onto their own periodic image. Everything else
   is the sol (unreacted monomer and unattached oligomers).
2. Dangling material: pieces of the gel attached through a single bond
   (a bridge) that lead to no infinite (winding) part, pruned repeatedly.
3. Active network: the rest. Junctions are its nodes with three or more
   active bonds; strands are the chains of two-bond nodes between them.
   A strand that returns to its own junction without winding is a loop and
   is not elastically effective.

Densities use the cell volume. M_c = rho N_A / nu_e. The rubbery shear
modulus estimates are the affine (nu_e k T) and phantom (xi k T / V,
xi = strands - junctions) limits; they ignore entanglements and are
for the rubbery state above Tg.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

NA = 6.02214076e23
KB = 1.380649e-23


def _lifts(nodes, edges, start):
    """Integer cell offsets of nodes reached from start, and whether a cycle winds.

    edges: node -> [(other, shift)] where shift is the cell offset of other
    seen from node across that edge.
    """
    lift, winding, pending = {start: np.zeros(3, dtype=int)}, False, [start]
    while pending:
        i = pending.pop()
        for j, shift in edges[i]:
            if j not in nodes:
                continue
            expected = lift[i] + shift
            if j not in lift:
                lift[j] = expected
                pending.append(j)
            elif (expected != lift[j]).any():
                winding = True
    return lift, winding


def _bridges(adjacency):
    """Edge ids that are bridges; adjacency: node -> [(other, edge id)]. Iterative Tarjan."""
    order, low, bridges, counter = {}, {}, set(), 0
    for root in adjacency:
        if root in order:
            continue
        order[root] = low[root] = counter; counter += 1
        stack = [(root, None, iter(adjacency[root]))]
        while stack:
            node, via, neighbors = stack[-1]
            advanced = False
            for other, edge in neighbors:
                if edge == via:
                    continue
                if other not in order:
                    order[other] = low[other] = counter; counter += 1
                    stack.append((other, edge, iter(adjacency[other])))
                    advanced = True
                    break
                low[node] = min(low[node], order[other])
            if not advanced:
                stack.pop()
                if stack:
                    parent = stack[-1][0]
                    low[parent] = min(low[parent], low[node])
                    if low[node] > order[parent]:
                        bridges.add(via)
    return bridges


def _backbone_bonds(atom_edges, start, end, interior):
    """Bonds on the shortest heavy-atom path from a start junction atom to an end junction atom.

    Only interior atoms of the strand are crossed. None when both ends are the same
    junction (a strand closing through the periodic boundary).
    """
    start, end = set(start), set(end)
    if start & end:
        return None
    frontier, seen, distance = list(start), set(start), 0
    while frontier:
        distance += 1
        following = []
        for a in frontier:
            for b, _ in atom_edges[a]:
                if b in end:
                    return distance
                if b in interior and b not in seen:
                    seen.add(b); following.append(b)
        frontier = following
    return None


def _ring_systems(graph, heavy, shift_of):
    """Union of small rings (3-8) that share atoms, over heavy atoms.

    A ring that closes only through the periodic boundary (its bond shifts do
    not sum to zero) is part of the infinite network, not a ring, and is skipped.
    """
    from polypaves.ring_piercing import RingFinder
    neighbors = [set() for _ in graph['atoms']]
    for i, j, _ in graph['bonds']:
        neighbors[i].add(j); neighbors[j].add(i)
    parent = list(range(len(graph['atoms'])))
    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for ring in RingFinder(neighbors).rings_near(heavy):
        if sum(shift_of(ring[m], ring[(m + 1) % len(ring)]) for m in range(len(ring))).any():
            continue
        for a in ring[1:]:
            parent[root(a)] = root(ring[0])
    return root


def polar_groups(graphs):
    """Counts of polar and reactive groups from the chemical graph."""
    counts = Counter()
    for graph in graphs:
        atoms = graph['atoms']
        neighbors = defaultdict(list)
        for i, j, order in graph['bonds']:
            neighbors[i].append((j, order)); neighbors[j].append((i, order))
        element = [a['atomic_number'] for a in atoms]
        from polypaves.ring_piercing import RingFinder
        rings = RingFinder({i: {j for j, _ in neighbors[i]} for i in range(len(atoms))}).rings_near(
            [i for i in range(len(atoms)) if element[i] != 1])
        epoxide_O = {a for r in rings if len(r) == 3 for a in r if element[a] == 8}
        aromatic_rings = sum(all(atoms[a]['aromatic'] for a in r) for r in rings)
        if aromatic_rings:
            counts['aromatic_rings'] += aromatic_rings
        for i, atom in enumerate(atoms):
            z = element[i]
            hydrogens = sum(element[j] == 1 for j, _ in neighbors[i])
            heavy = [(j, o) for j, o in neighbors[i] if element[j] != 1]
            if z == 8:
                if hydrogens >= 1:
                    aromatic = any(atoms[j]['aromatic'] for j, _ in heavy)
                    counts['phenolic_OH' if aromatic else 'aliphatic_OH'] += 1
                elif i in epoxide_O:
                    counts['epoxide_rings'] += 1
                elif len(heavy) == 2 and all(o == 1 for _, o in heavy):
                    counts['ether_O'] += 1
                elif len(heavy) == 1 and heavy[0][1] == 2:
                    counts['carbonyl_O'] += 1
            elif z == 7:
                counts[f'N_with_{hydrogens}H'] += 1
            if z in (7, 8):
                counts['hbond_acceptor_atoms'] += 1
                counts['hbond_donor_H'] += hydrogens
            if atom.get('reaction_state') == 'cationic_activity':
                counts['active_cationic_sites'] += 1
    return dict(counts)


def network_structure(graphs, positions_nm, box_nm, temperature_K=300., details=None):
    """Strand/junction statistics of a periodic network; positions and box in nm.

    details: optional list that receives one dict per strand or loop with the
    graph index and the atom indices of its two junctions and its interior.
    """
    xyz = np.asarray(positions_nm, dtype=float)
    box = np.asarray(box_nm, dtype=float)
    if box.shape == (3,):
        box = np.diag(box)
    inverse = np.linalg.inv(box)
    volume_m3 = abs(np.linalg.det(box)) * 1e-27
    masses = np.array([a['mass'] for g in graphs for a in g['atoms']], dtype=float)
    total_mass = masses.sum()

    sol_mass = dangling_mass = active_mass = 0.
    junction_degrees, strand_lengths, strand_masses, backbone_bonds = Counter(), [], [], []
    strands = loops = gel_molecules = 0
    offset = 0
    for gi, graph in enumerate(graphs):
        n = len(graph['atoms'])
        local_xyz, local_mass = xyz[offset:offset + n], masses[offset:offset + n]
        offset += n
        element = [a['atomic_number'] for a in graph['atoms']]
        heavy = [i for i in range(n) if element[i] != 1]
        # Hydrogen mass goes to its heavy atom.
        carried = local_mass.copy()
        bonds = []
        for i, j, _ in graph['bonds']:
            if element[i] == 1 or element[j] == 1:
                h, x = (i, j) if element[i] == 1 else (j, i)
                if element[x] != 1:
                    carried[x] += local_mass[h]; carried[h] = 0.
                continue
            bonds.append((i, j, -np.rint((local_xyz[j] - local_xyz[i]) @ inverse).astype(int)))
        if not heavy:
            sol_mass += local_mass.sum()
            continue
        # Contract ring systems; unwrap each system relative to its root.
        def shift_of(i, j):
            return -np.rint((local_xyz[j] - local_xyz[i]) @ inverse).astype(int)
        root = _ring_systems(graph, heavy, shift_of)
        atom_edges = defaultdict(list)
        for i, j, shift in bonds:
            atom_edges[i].append((j, shift)); atom_edges[j].append((i, -shift))
        members = defaultdict(list)
        for i in heavy:
            members[root(i)].append(i)
        offset_of = {}
        for node, atoms in members.items():
            lift, _ = _lifts(set(atoms), atom_edges, atoms[0])
            offset_of.update(lift)
        node_mass = {node: carried[atoms].sum() for node, atoms in members.items()}
        node_size = {node: len(atoms) for node, atoms in members.items()}
        edges, adjacency, edge_shift = defaultdict(list), defaultdict(list), {}
        for k, (i, j, shift) in enumerate(bonds):
            a, b = root(i), root(j)
            if a == b:
                continue
            s = shift + offset_of[i] - offset_of[j]
            edges[a].append((b, s)); edges[b].append((a, -s))
            adjacency[a].append((b, k)); adjacency[b].append((a, k))
            edge_shift[k] = (a, b, s)
        for node in members:
            adjacency[node]
        _, winds = _lifts(set(members), edges, root(heavy[0]))
        if not winds:
            sol_mass += local_mass.sum()
            continue
        gel_molecules += 1
        # 2-edge-connected components and whether each one winds.
        bridge_ids = _bridges(adjacency)
        component, comps = {}, []
        for node in members:
            if node in component:
                continue
            todo, comps_nodes = [node], []
            component[node] = len(comps)
            while todo:
                x = todo.pop(); comps_nodes.append(x)
                for y, k in adjacency[x]:
                    if k not in bridge_ids and y not in component:
                        component[y] = len(comps); todo.append(y)
            comps.append(comps_nodes)
        inner = defaultdict(list)
        for k, (a, b, s) in edge_shift.items():
            if k not in bridge_ids:
                inner[a].append((b, s)); inner[b].append((a, -s))
        winding = [_lifts(set(nodes), inner, nodes[0])[1] for nodes in comps]
        # Prune non-winding leaves of the bridge tree: dangling material.
        tree = defaultdict(set)
        for k in bridge_ids:
            a, b, _ = edge_shift[k]
            tree[component[a]].add((component[b], k)); tree[component[b]].add((component[a], k))
        alive = set(range(len(comps)))
        leaves = [c for c in alive if len(tree[c]) <= 1 and not winding[c]]
        while leaves:
            c = leaves.pop()
            if c not in alive:
                continue
            alive.discard(c)
            for other, k in list(tree[c]):
                tree[other].discard((c, k))
                if other in alive and len(tree[other]) <= 1 and not winding[other]:
                    leaves.append(other)
        active_nodes = {x for c in alive for x in comps[c]}
        dangling_mass += sum(node_mass[x] for x in members if x not in active_nodes)
        active_mass += sum(node_mass[x] for x in active_nodes)
        # Strands between junctions of the active graph.
        active = defaultdict(list)
        for k, (a, b, s) in edge_shift.items():
            if a in active_nodes and b in active_nodes:
                active[a].append((b, k, s)); active[b].append((a, k, -s))
        junctions = {x for x in active_nodes if len(active[x]) >= 3}
        for x in junctions:
            junction_degrees[len(active[x])] += 1
        used = set()
        for start in junctions:
            for other, k, s in active[start]:
                if k in used:
                    continue
                used.add(k)
                path_mass, path_atoms, shift, previous, node = 0., 0, s.copy(), k, other
                interior = []
                while node not in junctions:
                    path_mass += node_mass[node]; path_atoms += node_size[node]
                    interior.append(node)
                    (nxt, k2, s2), = [e for e in active[node] if e[1] != previous]
                    used.add(k2); shift += s2; previous, node = k2, nxt
                is_loop = node == start and not shift.any()
                if details is not None:
                    details.append(dict(graph=gi, loop=is_loop, junctions=(members[start], members[node]),
                                        interior=[members[x] for x in interior]))
                if is_loop:
                    loops += 1
                    continue
                strands += 1
                strand_lengths.append(path_atoms)
                strand_masses.append(path_mass)
                bonds_along = _backbone_bonds(atom_edges, members[start], members[node],
                                              {a for x in interior for a in members[x]})
                if bonds_along is not None:
                    backbone_bonds.append(bonds_along)
        # Junction-free active cycles (a chain wrapping the cell by itself).
        for x in active_nodes - junctions:
            for _, k, _ in active[x]:
                if k not in used:
                    ring_edges = {k}
                    todo = [x]
                    seen = {x}
                    while todo:
                        y = todo.pop()
                        for z, k3, _ in active[y]:
                            ring_edges.add(k3)
                            if z not in seen:
                                seen.add(z); todo.append(z)
                    used |= ring_edges
                    strands += 1
                    strand_lengths.append(sum(node_size[y] for y in seen))
                    strand_masses.append(sum(node_mass[y] for y in seen))
    junction_count = sum(junction_degrees.values())
    density_kg_m3 = total_mass / 1000 / NA / volume_m3
    nu_e = strands / NA / volume_m3            # mol/m^3
    xi = (strands - junction_count) / NA / volume_m3
    lengths = np.asarray(strand_lengths, dtype=float)
    polar = polar_groups(graphs)
    per_volume = 1 / NA / volume_m3 / 1000       # count -> mmol/cm^3 (= mol/L)
    return dict(
        definitions='heavy-atom graph, fused 3-8 rings contracted; gel = periodically winding molecules; '
                    'dangling = pieces behind a bridge with no winding part; junction = active node with >= 3 active bonds',
        volume_nm3=volume_m3 * 1e27, density_g_cm3=density_kg_m3 / 1000,
        gel_molecules=gel_molecules,
        sol_mass_fraction=float(sol_mass / total_mass),
        dangling_mass_fraction=float(dangling_mass / total_mass),
        active_mass_fraction=float(active_mass / total_mass),
        junctions=junction_count,
        junction_functionality=dict(sorted((int(k), v) for k, v in junction_degrees.items())),
        mean_junction_functionality=float(sum(k * v for k, v in junction_degrees.items()) / junction_count)
            if junction_count else None,
        effective_strands=strands, elastically_ineffective_loops=loops,
        strand_density_nu_e_mol_m3=float(nu_e),
        strand_density_nu_e_mmol_cm3=float(nu_e / 1000),
        junction_density_mmol_cm3=float(junction_count * per_volume),
        cycle_rank_density_mmol_cm3=float(xi / 1000),
        Mc_g_mol=float(density_kg_m3 * 1000 / nu_e) if strands else None,
        strand_heavy_atoms=dict(mean=float(lengths.mean()), median=float(np.median(lengths)),
                                min=int(lengths.min()), max=int(lengths.max())) if strands else None,
        mean_strand_interior_mass_g_mol=float(np.mean(strand_masses)) if strands else None,
        # Every effective strand: interior heavy atoms (ring systems count all their atoms) -> count.
        strand_heavy_atoms_histogram={int(k): v for k, v in sorted(Counter(strand_lengths).items())},
        # Contour: bonds on the shortest path between the two junction atoms (side atoms of rings excluded).
        strand_backbone_bonds=dict(mean=float(np.mean(backbone_bonds)), median=float(np.median(backbone_bonds)),
                                   min=int(min(backbone_bonds)), max=int(max(backbone_bonds))) if backbone_bonds else None,
        strand_backbone_bonds_histogram={int(k): v for k, v in sorted(Counter(backbone_bonds).items())},
        rubbery_shear_modulus_estimates_MPa=dict(
            temperature_K=temperature_K,
            affine=float(nu_e * KB * NA * temperature_K / 1e6),
            phantom=float(xi * KB * NA * temperature_K / 1e6)),
        polar_groups_count=polar,
        polar_groups_mmol_cm3={k: float(v * per_volume) for k, v in polar.items()},
        polar_groups_mmol_g={k: float(v / total_mass * 1000) for k, v in polar.items()},
        limitations=['Topological counts; entanglements and trapped entanglements are not included.',
                     'Rubbery modulus estimates apply above Tg, not to the glass.',
                     'Ring systems up to 8-membered are contracted; larger cycles are network loops.'])


def analyze_file(path, temperature_K=300.):
    """Analyze a saved *.structure.json (positions and box in angstrom)."""
    data = json.loads(Path(path).read_text())
    graphs = data['graphs']
    host = data.get('host_atoms')
    if host:  # a wet snapshot: keep the network only
        kept, count = [], 0
        for g in graphs:
            if count >= host:
                break
            kept.append(g); count += len(g['atoms'])
        graphs = kept
    n = sum(len(g['atoms']) for g in graphs)
    return network_structure(graphs, np.asarray(data['positions'])[:n] / 10.,
                             np.asarray(data['box_vectors']) / 10., temperature_K)


def plot_strand_histogram(result, path, title=''):
    """Strand length histograms: backbone bonds (contour) and interior heavy atoms (size)."""
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
    for ax, key, label in ((axes[0], 'strand_backbone_bonds_histogram', 'Backbone bonds between junctions'),
                           (axes[1], 'strand_heavy_atoms_histogram', 'Interior heavy atoms (whole ring systems)')):
        histogram = {int(k): v for k, v in result[key].items()}
        ax.bar(list(histogram), list(histogram.values()), width=.8)
        ax.set(xlabel=label, ylabel='Strands')
    fig.suptitle(title or f"{result['effective_strands']} effective strands, Mc {result['Mc_g_mol']:.0f} g/mol")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_strand_csv(result, path):
    """Both histograms in one table: length, strands by backbone bonds, strands by interior heavy atoms."""
    bonds = {int(k): v for k, v in result['strand_backbone_bonds_histogram'].items()}
    atoms = {int(k): v for k, v in result['strand_heavy_atoms_histogram'].items()}
    with Path(path).open('w') as out:
        out.write('length,strands_by_backbone_bonds,strands_by_interior_heavy_atoms\n')
        for length in sorted(set(bonds) | set(atoms)):
            out.write(f'{length},{bonds.get(length, 0)},{atoms.get(length, 0)}\n')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Network graph analysis of a saved *.structure.json')
    parser.add_argument('structure', type=Path)
    parser.add_argument('--histogram-csv', type=Path, help='write the strand length histogram as CSV')
    parser.add_argument('--histogram-png', type=Path, help='plot the strand length histogram')
    args = parser.parse_args()
    result = analyze_file(args.structure)
    if args.histogram_csv:
        write_strand_csv(result, args.histogram_csv)
    if args.histogram_png and result['effective_strands']:
        plot_strand_histogram(result, args.histogram_png)
    print(json.dumps(result, indent=2))
