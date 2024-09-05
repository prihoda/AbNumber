#!/usr/bin/env python

import argparse
import os
import sys
import pandas
from abnumber import Position

def merge_anarci_csv(input_paths, output_path, scheme, add_regions=False):
    assert input_paths, 'No input files specified'

    # Merge position numbers
    metadata_columns = list()
    position_numbers = set()
    is_heavy_chain = None
    for path in input_paths:
        df = pandas.read_csv(path, nrows=1)
        assert 'chain_type' in df.columns, f'Column chain_type not found in {path}'
        if is_heavy_chain is None:
            is_heavy_chain = df['chain_type'].iloc[0] == 'H'
        else:
            assert is_heavy_chain == (df['chain_type'].iloc[0] == 'H'), f'All files need to have the same chain type, inconsistent chain type in {path}'
        df_positions = df.loc[:, ('1' if '1' in df.columns else '2'):].columns.tolist()
        position_numbers.update(df_positions)
        for c in df.columns:
            if c not in position_numbers and c not in metadata_columns:
                metadata_columns.append(c)
    # Use the correct sorting by using Position objects
    sorted_positions = sorted([Position.from_string(p, 'H' if is_heavy_chain else 'L', scheme) for p in position_numbers])
    sorted_position_numbers = [p.format(chain_type=False) for p in sorted_positions]

    # Verify that all positions are present
    missing_positions = set(position_numbers) - set(sorted_position_numbers)
    if missing_positions:
        raise ValueError(f'Unexpected error merging positions, missing: {", ".join(missing_positions)}')

    # Merge files
    num = 0
    with open(output_path, 'w') as f:
        # Write header
        if add_regions:
            f.write(','.join(metadata_columns + [p.get_region()+'_'+n for p, n in zip(sorted_positions, sorted_position_numbers)]) + os.linesep)
        else:
            f.write(','.join(metadata_columns + sorted_position_numbers) + os.linesep)

        # Write files one by one
        for i, path in enumerate(input_paths):
            df = pandas.read_csv(path)
            df_positions = df.loc[:, ('1' if '1' in df.columns else '2'):].columns.tolist()
            assert df_positions == [c for c in sorted_position_numbers if c in df_positions], \
                f'Unexpected reordering of positions in CSV file: {sorted_position_numbers}'

            df = df.reindex(metadata_columns + sorted_position_numbers, axis=1)
            df.loc[:, sorted_position_numbers].fillna('-', inplace=True)
            df.to_csv(f, index=False, header=False)
            num += len(df)

    return num

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge ANARCI CSV files')
    parser.add_argument('--scheme', help='Numbering scheme')
    parser.add_argument('--add-regions', help='Add region prefixes to header', action='store_true')
    parser.add_argument('input_paths', nargs='+', help='Paths to input CSV files')
    parser.add_argument('output_path', help='Path to output CSV file')
    args = parser.parse_args()

    if os.path.exists(args.output_path):
        raise ValueError(f'Output file {args.output_path} already exists')

    print(f'Merging {len(args.input_paths):,} files...', file=sys.stderr)
    try:
        num = merge_anarci_csv(args.input_paths, args.output_path, args.scheme, add_regions=args.add_regions)
        print(f'Merged {num:,} sequences from {len(args.input_paths):,} files', file=sys.stderr)
        print(f'Saved output to: {args.output_path}', file=sys.stderr)
    except Exception as e:
        # Remove output file if an error occurred
        if os.path.exists(args.output_path):
            os.remove(args.output_path)
        raise
