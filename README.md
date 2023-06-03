# Nextflow ONT Influenza Typing Pipeline
A nextflow pipeline built to type influenza sequencing data. The pipeline was inspired by the [INSaFLU Database](https://insaflu.insa.pt/) and the [wf-flu pipeline](https://github.com/epi2me-labs/wf-flu) from EPI2ME Labs. Both of these tools require a connection to the internet, which is not possible in many cases. Thus, this pipeline is a quick alternative that can be run locally.

The pipeline processes/filters reads, aligns reads to the INSaFLU database using [Minimap2](https://github.com/lh3/minimap2), and selects the gene variants with the highest number of reads to assemble to (you can also set the pipeline to assemble to any gene variant above a certain read threshold - See **Additional Options**). Then, it performs a **reference-based** assembly to these sequences using Minimap2 for alignment, [Medaka](https://github.com/nanoporetech/medaka) and [Longshot](https://github.com/pjedge/longshot) for variant calling, and [BCFTools](https://samtools.github.io/bcftools/) and [BEDTools](https://github.com/arq5x/bedtools2) to generate consensus sequences. Finally, the assembled sequences are typed by aligning them against the INSaFLU typing database using [Abricate](https://github.com/tseemann/abricate).

## Installation

To install the pipeline, enter the following commands:
```
# Clone the repository
git clone https://github.com/rchapman2000/ont-influenza-typing

# Create a conda environment using the provided environment.yml file
conda env create -f environment.yml

# Activate the conda environment
conda activate ONT-Flu-Typing
```

### Updating the Pipeline
If you already have the pipeline installed, you can update it using the following commands:
```
# Navigate to your installation directory
cd ont-influenza-typing

# Use git pull to get the latest update
git pull

# Activate the conda environment and use the environment.yml file to download updates
conda activate ONT-Flu-Typing
conda env update --file environment.yml --prune
```

## Usage
To run the pipeline, use the following command:
```
# You must either be in the same directory as the main.nf file or reference the file location.
nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR --model MEDAKA_MODEL
```

### Options
The pipeline also supports the following optional arguments:

| Option | Type | Description |
|---|---|---|
| --trimONTAdapters | *None* | Enables ONT Adapter/Barcodee trimming using Porechop [Default = OFF] |
| --minReadLen | *INT* | If supplied, the pipeline will perform length filtering using Chopper excluding reads less than this size and sets the minimum length for Miniasm assembly [Default = 200 bp] |
| --maxReadLen | *INT* | If supplied, the pipeline will perform legnth filtering using Chopper excluding reads greater than this size [Default = off]  |
| --minReadsToAssemble | *INT* | By default, the pipeline will assemble to the HA/NA segments with the most reads aligned to them. Supplying this parameter changes this to assemble reads to any segment with more than the supplied number of reads [Default = OFF]|
| --minCov | *INT* | The minimum coverage below which a position will be masked [Default = 10] |

To view the list of options from the command line, use the following command:
```
nextflow run main.nf --help
```
