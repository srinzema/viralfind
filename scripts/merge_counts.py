import sys
import pandas as pd
from pathlib import Path

out_file = sys.argv[1]
print(out_file)
count_table = pd.DataFrame()
for file in sys.argv[2:]:
    file = Path(file)
    sample = file.name.replace(".tsv", "")

    _df = pd.read_csv(file, sep="\t", header=1, index_col=0)
    _df = _df.drop(labels=["Chr", "Start", "End", "Strand", "Length"], axis=1)
    _df = _df.rename({_df.columns[0]: sample}, axis=1)

    count_table = pd.concat([count_table, _df], axis=1)

count_table.to_csv(out_file, sep="\t")
