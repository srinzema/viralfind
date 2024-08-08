import genomepy
import pandas as pd
import argparse
import time
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from loguru import logger
from pathlib import Path

FILE_EXTENSIONS = (
        "annotation.bed",
        "annotation.gtf",
        "fa",
        "fa.fai",
        "fa.sizes",
        "gaps.bed"
    )

def main(assembly_file, outdir, outname, cores, tmp=None):
    if not os.path.isfile(assembly_file):
        raise FileNotFoundError(f"The file {assembly_file} does not exist.")
    
    assemblies = load_assemblies(assembly_file)
    print(assemblies.head())

    # Install all the genomes
    assembly_directories = install_assemblies(assemblies, tmp or outdir, cores)
    create_metagenome(assembly_directories, outdir, outname)

    # Create the new directory in outdir
    new_dir_path = os.path.join(outdir, outname)
    if not os.path.isdir(new_dir_path):
        os.makedirs(new_dir_path)
    

def load_assemblies(assembly_file: str) -> pd.DataFrame:
    _df = pd.read_table(assembly_file)
    original_size = len(_df)
    # Drop all rows where annotation is False
    _df = _df[_df.annotation]
    print(f"Dropped {original_size - len(_df)} rows without annotation.")

    for column in _df.columns:
        _df[column] = _df[column].apply(lambda x: str(x).replace(" ", "_"))
    return _df


def install_assemblies(assemblies, working_dir, cores, wait_time=1):
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    print(f"Using: {working_dir}")

    genome_directories = []
    for _, row in assemblies.iterrows():
        name = row["name"]
        provider = row["provider"]
        localname = row["species"]
        assembly_directory = Path(working_dir) / localname
        genome_directories.append(assembly_directory)

        missing_files = False
        for extension in FILE_EXTENSIONS:
            file: Path = assembly_directory / f"{localname}.{extension}"
            if not file.exists():
                missing_files = True
                break

        if not missing_files:
            continue
        
        print(f"Installing: {name} at {assembly_directory}")
        try:
            # genomepy.clean()  # Clean before each install because of stupid "provider is offline" exceptions that I haven't been able to catch yet because of the weird design of genomepy.
            genomepy.install_genome(
                name = name,
                provider = provider,
                localname = localname,
                genomes_dir = working_dir,
                annotation = True,
                threads = cores
            )
        except Exception as e:
            print(e)

        time.sleep(wait_time)

    print("Downloaded all genomes.")
    print(f"Paths used: {genome_directories}")
    return sorted(genome_directories)


def create_metagenome(assembly_directories, outdir, genome_name):
    outdir = Path(outdir) / genome_name
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Making metagenome at: {outdir}")

    for extension in FILE_EXTENSIONS:
        out_file = outdir / f"{genome_name}.{extension}"
        open(out_file, "w").close()

        print(f"Constructing {out_file}")
        for assembly in assembly_directories:
            in_file = assembly / f"{assembly.name}.{extension}"
            append_file_contents(in_file, out_file, assembly.name)


def append_file_contents(input_file, output_file, origin_genome):
    """
    Appends the contents of input_file to output_file.
    
    Args:
        input_file (str): Path to the file to read from.
        output_file (str): Path to the file to append to.
    """
    try:
        with open(input_file, 'r') as infile:
            with open(output_file, 'a') as outfile:
                outfile.write(
                    infile.read().replace(
                        'gene_id "', f'gene_id "{origin_genome}:'
                    )
                )
        print(f"Appended contents of {input_file} to {output_file}")
    except FileNotFoundError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some file paths and create directories.")
    parser.add_argument('assembly_file', type=str, help="Path to the assembly file (TSV format).")
    parser.add_argument('outdir', type=str, help="Path to the output directory where the new directory will be placed.")
    parser.add_argument('outname', type=str, help="Name of the new directory to be created within outdir.")
    parser.add_argument('--tmp', type=str, help="Path to the temporary directory. If not provided, outdir is used.", default=None)
    parser.add_argument('--cores', type=int, help="Number of CPU cores to use.", default=8)


    args = parser.parse_args()
    main(args.assembly_file, args.outdir, args.outname, args.cores, args.tmp)