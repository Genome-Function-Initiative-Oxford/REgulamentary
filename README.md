# REgulamentary: de novo genome-wide annotation of cis-regulatory elements
### Enhancer, Promoter, and CTCF

***
## Getting started
REgulamentary can be run using the __re__ conda environment, please follow the installation instruction detailed below. In doing so, the anlaysis is highly reproducible. 

### Installation instructions for conda environment

#### 1. Clone the repository
```
git clone git@github.com:Genome-Function-Initiative-Oxford/REgulamentary.git
cd REgulamentary
```

#### 2. Install anaconda
Check if [Anaconda](https://www.anaconda.com), [Miniconda](https://docs.conda.io/en/latest/miniconda.html), or [Mambaforge](https://mamba.readthedocs.io/en/latest/installation.html) is installed, using:
```
which conda
```   
If installed, the output should be:
```
~/anaconda3/condabin/conda
```
If [Anaconda](https://www.anaconda.com), [Miniconda](https://docs.conda.io/en/latest/miniconda.html), or [Mambaforge](https://mamba.readthedocs.io/en/latest/installation.html) is not installed, we recommend to install [Mambaforge](https://mamba.readthedocs.io/en/latest/installation.html), since it has already integrated mamba for a fast and parallelisable installation.   


Download [Mambaforge](https://mamba.readthedocs.io/en/latest/installation.html) ([Anaconda](https://www.anaconda.com) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)):
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```
- Run the installer as follows, and follow the on-screen commands.
```
sh Mambaforge-Linux-x86_64.sh
``` 
#### 3. Create Anaconda environment
Activate the conda 'base' environment (if not active): 
```
conda activate base
```

There are two ways to create the conda env upstream environment:
1) Using mamba (if [Mambaforge](https://mamba.readthedocs.io/en/latest/installation.html) was installed), and follow the on screen instructions:
```
mamba env create --file=envs/re.yml
```
2) Using conda, and follow the on screen instructions.
```
conda env create --file=envs/re.yml
```

#### 4. Activate the environment
Now, the upstream environment is created it needs to be activated: 
```
conda activate re
```
You can then use all of our upstream pipelines using this environment, enjoy!

### Environment installation note
REgulamentary has been successfully tested for the following operating systems: Ubuntu, CentOS, macOS (Intel CPU), and Windows. Unfortunately, it is not possible to install on macOS with M CPUs at the moment. 
For any error in the installation step, please open an [issue](https://github.com/Genome-Function-Initiative-Oxford/REgulamentary/issues) so we can give a general solution for users.

### Reproducibility :repeat:
If required for publication, package versions within the environment can be exported as follows:
```
conda env export > re_environment_versions.yml
```

***

### Pipeline updates :construction:
If any changes are made to REgulamentary workflow, it is possible to update the repository by entering the main folder and pulling the update using:
   ```
   # Enter the main folder
   cd REgulamentary

   # Pull updates
   git pull           
   ```
Alternatively, remove the cloned repository and then re-clone the repository as described above.   
Warning: use rm carefully!

```
rm -rf REgulamentary
``` 
<hr>

### :warning: Warning for University of Oxford CCB users :warning:
When using this repository, use the default terminal and __do not__ load any module in the server (if logged-in).

***

### Contact us
If you have any suggestions, spot any errors, or have any questions regarding REgulamentary, please do no hesitate to contact us anytime.   

:email: &emsp; [<simone.riva@imm.ox.ac.uk>](simone.riva@imm.ox.ac.uk)


### Running with BigWig Files

The pipeline can be run using BigWig files instead of BAMs by utilizing the signal value at each region instead of read counts. While results may differ slightly from those obtained using BAM files, if the BigWig files were generated at the highest possible resolution, the results are typically nearly identical.

An example configuration for running in this mode is provided at `config/astrocyte_chr21_bw.yaml`. To use BigWigs, simply specify the file paths for each mark (CTCF, H3K4me1, H3K4me3, H3K27ac) and ATAC/DNase under the `bigwigs` parameter. 

If you want the pipeline to automatically call peaks, add `create_bed_files: True` under the `bigwigs` section. Otherwise, you must provide BED files manually via the `union_peaks` parameter.

**Example:**

```yaml
bigwigs:
  ATAC: "test_data/bigwigs/astrocyte_DNase_chr21.bw"
  H3K4me1: "test_data/bigwigs/astrocyte_H3K4me1_chr21.bw"
  H3K4me3: "test_data/bigwigs/astrocyte_H3K4me3_chr21.bw"
  H3K27ac: "test_data/bigwigs/astrocyte_H3K27ac_chr21.bw"
  CTCF: "test_data/bigwigs/astrocyte_CTCF_chr21.bw"
  create_bed_files: True
```

To test this mode, ensure you are in the main project directory and run:

```bash
snakemake --configfile=config/astrocyte_chr21_bw.yaml all --cores 2