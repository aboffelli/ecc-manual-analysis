# Workflow of the TMA scripts

## QuPath/Stardist nuclear detection
Follow the instructions to install [QuPath](https://qupath.github.io/) and 
the [Stardist extension](https://github.com/qupath/qupath-extension-stardist).  
We use a pre-trained model available for download in [here](https://github.com/qupath/models/tree/main/stardist).

### Running the detection

The first step is to create a project for the TMA, and import the whole TMA image in QuPath.

Open the script ***01_tma_dearrayer.groovy*** in the Script editor (Automate > Script editor), and click Run. This script will detect the TMA columns and rows and save each core as a separate image.

Import all the core images to the QuPath project.

Open the script ***02_stardist_script.groovy*** and click **Run for project** (on the three dots beside the Run button). Select all your images. Remember to remove the full TMA image from the list, if you did not remove from the project.

When the script finished running for all images, open the script ***03_dab_export_measurements.groovy***, and click Run. This will generate a table containing all the information of every detection in each core.

Repeat the process for each TMA.

## Snakemake workflow

### Installing Snakemake
Ensure you have Snakemake installed. You can install it via conda or mamba:

```sh
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

Alternatively, install it using pip:
```sh
pip install snakemake
```

### Configure the workflow

Edit the ```config.yaml``` file to specify the input data and output directory.

```yaml
input_files:
  manual_scores: "path/to/data.txt"
  patient_file: "path/to/data.txt"
  tma_map: "path/to/data.txt"
  immuno_tables: "path/to/data.txt"
  original_gex: "path/to/data.txt"
  detection_tables: "path/to/data.txt"

output_dir: results
```

### Run snakemake

To test if everything is in place, we can do a dry run
```sh
snakemake -n
```

If the dry run returns no errors, we can execute the workflow
```sh
# Adjust the number of cores according to your preferences
snakemake --cores 6
```
