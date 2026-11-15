"""Kinetic chain length per active end, traced through the activity transfers.

Activation puts the activity state on one atom per start (the OH oxygen of
the opened ring). Every growth reaction moves it from one atom (old state
set, new state cleared) to another (old cleared, new set): the OH* of the
chain opens the next ring and the activity passes to that ring's oxygen.
Following that transfer assigns each reaction to the start it came from.

Chain length = epoxide units per active end = 1 (the ring opened by the
activation) + the propagation reactions of that end.
"""
import csv
import json
from pathlib import Path

import numpy as np


def _transfers(changes, state):
    lost = [sid for sid, c in changes.items()
            if c.get('old', {}).get('reaction_state') == state and c.get('new', {}).get('reaction_state') != state]
    gained = [sid for sid, c in changes.items()
              if c.get('new', {}).get('reaction_state') == state and c.get('old', {}).get('reaction_state') != state]
    return lost, gained


def kinetic_chains(activation_records, reaction_records, state='cationic_activity'):
    """Per-start chain records from the activation and reaction histories (in order)."""
    owner, chains = {}, []
    for record in activation_records:
        _, gained = _transfers(record['reaction']['atom_changes'], state)
        for sid in sorted(gained):
            owner[sid] = len(chains)
            chains.append(dict(chain=len(chains), start_atom=sid, propagations=0, active_atom=sid))
    for record in reaction_records:
        lost, gained = _transfers(record['atom_changes'], state)
        if len(lost) != 1 or len(gained) != 1 or lost[0] not in owner:
            raise ValueError(f'reaction in cycle {record.get("cycle")} does not move one tracked activity '
                             f'(lost {lost}, gained {gained})')
        chain = owner.pop(lost[0])
        owner[gained[0]] = chain
        chains[chain]['propagations'] += 1
        chains[chain]['active_atom'] = gained[0]
    for chain in chains:
        chain['epoxide_units'] = chain['propagations'] + 1
    return chains


def chain_statistics(chains):
    units = np.array([c['epoxide_units'] for c in chains], dtype=float)
    if not len(units):
        return dict(active_ends=0)
    values, counts = np.unique(units.astype(int), return_counts=True)
    number = float(units.mean())
    weight = float((units ** 2).sum() / units.sum())
    return dict(active_ends=len(units), propagation_reactions=int(units.sum() - len(units)),
                number_average_epoxide_units=number, weight_average_epoxide_units=weight,
                dispersity=weight / number, std_epoxide_units=float(units.std()),
                min_epoxide_units=int(units.min()), median_epoxide_units=float(np.median(units)),
                max_epoxide_units=int(units.max()),
                distribution=[dict(epoxide_units=int(v), active_ends=int(n)) for v, n in zip(values, counts)],
                definition='epoxide units per active end = 1 (ring opened by the activation) + growth reactions '
                           'traced through the activity transfer; number/weight averages over active ends')


def read_jsonl(path):
    path = Path(path)
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def write_chain_lengths(output, state='cationic_activity'):
    """Chain lengths of a finished crosslinking run: JSON summary, per-chain CSV and a histogram PNG."""
    output = Path(output)
    chains = kinetic_chains(read_jsonl(output / 'activation_history.jsonl'),
                            read_jsonl(output / 'reaction_history.jsonl'), state)
    summary = chain_statistics(chains)
    (output / 'chain_lengths.json').write_text(json.dumps(summary, indent=2) + '\n')
    with open(output / 'chain_lengths.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['chain', 'start_atom', 'active_atom', 'propagations', 'epoxide_units'])
        writer.writeheader()
        writer.writerows(chains)
    if chains:
        try:
            import matplotlib
            matplotlib.use('Agg')
            from matplotlib import pyplot as plt
            fig, ax = plt.subplots(figsize=(5.5, 3.6))
            values = [d['epoxide_units'] for d in summary['distribution']]
            ax.bar(values, [d['active_ends'] for d in summary['distribution']], width=.8, color='tab:blue')
            ax.axvline(summary['number_average_epoxide_units'], c='tab:red', ls='--', lw=1,
                       label=f"number average {summary['number_average_epoxide_units']:.2f}")
            ax.axvline(summary['weight_average_epoxide_units'], c='tab:orange', ls=':', lw=1.2,
                       label=f"weight average {summary['weight_average_epoxide_units']:.2f}")
            ax.set(xlabel='Epoxide units per active end', ylabel='Active ends',
                   title=f"Kinetic chain length ({summary['active_ends']} active ends)")
            ax.legend(fontsize=8)
            fig.tight_layout()
            fig.savefig(output / 'chain_lengths.png', dpi=160)
            plt.close(fig)
        except Exception as exc:  # a figure must never fail a finished build
            print(f'WARNING: chain length figure not written: {exc}', flush=True)
    return summary
