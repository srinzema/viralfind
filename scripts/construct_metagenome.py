import genomepy
import pandas as pd
import argparse
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

def main(assembly_file, outdir, outname, cores, tmp=None):
    # Check if assembly_file exists
    if not os.path.isfile(assembly_file):
        raise FileNotFoundError(f"The file {assembly_file} does not exist.")
    
    assemblies = load_assemblies(assembly_file)
    print(assemblies.head())

    # Install all the genomes
    assembly_directories = install_assemblies(assemblies, tmp or outdir, cores)
    print(f"Assemblies: {assembly_directories}")
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


def install_assemblies(assemblies, working_dir, cores):
    def install_genome_worker(name, provider, localname, genomes_dir):
        print(f"Installing: {name} at {genomes_dir}/{localname}")
        try:
            genomepy.install_genome(
                name=name,
                provider=provider,
                localname=localname,
                genomes_dir=genomes_dir,
                annotation=True,
                force=True,
            )
        except Exception as e:
            print(f"Exception: {e}")

    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    print(f"Using: {working_dir}")

    futures = []
    genome_directories = []
    # Use ThreadPoolExecutor to handle multithreading
    with ThreadPoolExecutor(max_workers=cores) as executor:
        for _, row in assemblies.iterrows():
            # Submit tasks to the thread pool
            future = executor.submit(
                install_genome_worker,
                row["name"],
                row["provider"],
                row["species"],
                working_dir
            )
            futures.append(future)
            genome_directories.append(
                Path(working_dir) / row["species"]
            )

        # Wait for all futures to complete
        for future in as_completed(futures):
            # Optionally handle exceptions or results here
            future.result()

    print("Downloaded all genomes.")
    return sorted(genome_directories)


def create_metagenome(assembly_directories, outdir, genome_name):
    outdir = Path(outdir) / genome_name
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Making metagenome at: {outdir}")

    file_extensions = (
        "annotation.bed",
        "annotation.gtf",
        "fa",
        "fa.fai",
        "fa.sizes",
        "gaps.bed"
    )

    for extension in file_extensions:
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