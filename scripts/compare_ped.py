#!/usr/bin/env python3

from cloudpathlib import AnyPath
from collections import defaultdict
import click
import json


@click.command()
@click.option(
    '--king_relatedness',
    help='The cuKING relatedness JSON output file path.',
    required=True,
)
@click.option(
    '--king_threshold',
    default=0.177,  # 2nd degree
    help='The coefficient value cut-off to use to only filter closely related samples.',
)
@click.option(
    '--pedigree', help='A pedigree file name to compare against.', required=True
)
def main(king_relatedness, king_threshold, pedigree):
    """Compares a cuKING relatedness output with a given pedigree."""
    king_related = defaultdict(set)
    with AnyPath(king_relatedness).open('r') as f:
        king_map = json.load(f)
        for id, entries in king_map.items():
            for related, coeff in entries.items():
                if coeff < king_threshold:
                    continue
                king_related[id].add(related)
                king_related[related].add(id)

    ped_related = defaultdict(set)
    all_ids = []
    with AnyPath(pedigree).open('r') as f:
        for line in f.readlines()[1:]:  # Skip header
            entries = line.split()
            if len(entries) != 4:
                continue
            all_ids.append(entries[0])
            for parent in entries[1:3]:
                if parent != '0':
                    ped_related[entries[0]].add(parent)
                    ped_related[parent].add(entries[0])

    for id in all_ids:
        missing = ped_related[id] - king_related[id]
        if missing:
            print(f'Missing relatedness for {id}: {missing}')
        unexpected = king_related[id] - ped_related[id]
        if unexpected:
            print(f'Unexpected relatedness for {id}: {unexpected}')


if __name__ == '__main__':
    main()
