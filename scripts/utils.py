from pathlib import Path
import pandas as pd
import os, sys


def load_assemblies(assembly_file: str) -> pd.DataFrame:
    assemblies = pd.read_table(assembly_file)
    original_size = len(assemblies)
    # Drop all rows where annotation is False
    assemblies = assemblies[assemblies.annotation]
    print(f"Dropped {original_size - len(assemblies)} rows without annotation.")

    for column in assemblies.columns:
        assemblies[column] = assemblies[column].apply(
            lambda x: str(x).replace(" ", "_")
        )
    return assemblies


def load_samples(sample_file: str, fastq_dir: str) -> pd.DataFrame:
    samples = pd.read_table(sample_file, comment="#").set_index("samples", drop=True)
    if "alias" not in samples.columns:
        samples["alias"] = samples.index

    read1 = []
    read2 = []
    fastq_dir = Path(fastq_dir)
    for sample in samples.index:
        fastqs = [str(x) for x in fastq_dir.glob(f"{sample}*")]
        if len(fastqs) == 1:
            # one file found
            read1.append(fastqs[0])
            read2.append(pd.NA)

        elif len(fastqs) == 2:
            # two files found
            r1 = next(filter(lambda x: "R1" in x, fastqs))
            read1.append(r1)

            r2 = next(filter(lambda x: "R2" in x, fastqs))
            read2.append(r2)
        else:
            raise ValueError(
                f"Incorrect number of fastq files found for sample {sample} ({len(fastqs)})"
            )

    samples.insert(0, "read1", read1)
    samples.insert(1, "read2", read2)
    return samples
