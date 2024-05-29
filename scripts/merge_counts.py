import sys
import argparse
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(prog="merge_counts")
parser.add_argument("-o", "--out")
parser.add_argument("-s", "--samplesheet")
parser.add_argument("-m", "--merge", action="store_true")
parser.add_argument("-f", "--filter", action="store_true")
parser.add_argument("-i", "--input", nargs="+", default=[])
args = parser.parse_args()


samplesheet = pd.read_table(args.samplesheet).set_index("samples", drop=True)
if "alias" not in samplesheet.columns:
    samplesheet["alias"] = samplesheet.index
out_file = args.out
print(f" saving to {out_file}")

# Construct initial dataframe
files = args.input
count_table = pd.DataFrame()
for file in files:
    file = Path(file)
    sample = file.name.replace(".tsv", "")

    _df = pd.read_csv(file, sep="\t", header=1, index_col=0)
    _df = _df.drop(labels=["Chr", "Start", "End", "Strand", "Length"], axis=1)
    _df = _df.rename({_df.columns[0]: sample}, axis=1)

    count_table = pd.concat([count_table, _df], axis=1)

# Merge replicate groups
if args.merge and "replicate_group" in samplesheet.columns:
    for group in samplesheet["replicate_group"].unique():
        # Get columns for all samples in replicate group
        sample_names = samplesheet.loc[
            samplesheet["replicate_group"] == group
        ].alias.values
        sample_counts = count_table[sample_names]

        # Add a new column for the group containing the sum
        count_table[group] = sample_counts.sum(axis=1)

        # Drop old columns
        count_table = count_table.drop(sample_names, axis=1)

# Drop rows where sum is 0
if args.filter:
    count_table = count_table[count_table.sum(axis=1) > 0]

print(count_table.head())
count_table.to_csv(out_file, sep="\t")
