# 📦 Installation

## System requirements

**Operating system**

* Linux (x86_64, preferred and used for development and testing)
* MacOS (Intel or Apple Silicon; not rigorously tested, some bioconda dependencies might not work out-of-the-box, see [Troubleshooting](../principles/troubleshooting.md#working-on-apple-silicon))

**Software**

* [git](https://git-scm.com/)
* Conda, e.g., via [miniforge](https://github.com/conda-forge/miniforge) (recommended) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html)
* All other dependencies (Snakemake, Python and R packages) are installed through the conda environments under `envs/`, which list the dependencies and their version constraints.

**Hardware**

* The [demo](quickstart.md) runs on a standard desktop or laptop computer (CPU only, 8 GB RAM is sufficient).
* Disk space: ~8 GB for the conda environments of the demo (~5.5 GB environments, ~2.5 GB conda package cache).
* **Optional, non-standard hardware:** an NVIDIA GPU with CUDA 12 support.
  A GPU is only needed for deep-learning-based methods (e.g. scVI, scANVI, scPoli, DRVI, sysVI via the `scvi-tools` and `scarches` environments) and for GPU-accelerated preprocessing and clustering (`rapids_singlecell` environment, used when `use_gpu: true`).
  All other modules run on CPU (see [Working with CPUs only](../principles/troubleshooting.md#working-with-cpus-only)).
* For large atlases (millions of cells), we recommend a compute cluster; see {ref}`cluster_execution`.

## Clone the repository

Depending on whether you have set up SSH or HTTPS with [PAT](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens), you can clone the repository.

SSH:
```
git clone git@github.com:HCA-integration/scAtlasTb.git
```

HTTPS:
```
git clone https://github.com/HCA-integration/scAtlasTb.git
```

All following commands are run from the root of the cloned repository.

```
cd scAtlasTb
```

## Install dependencies

The modules are tested and developed using task-specific conda environments, which should be quick to set up when using [libmamba](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community).

> 📝  **Note** If you use conda version 22.11 or above, make sure you set the conda solver to [`libmamba`](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) for significantly faster installation. For newer versions or if you are using mamba directly, `libmamba` should already be the default.

All the conda environments used by the toolbox are under `envs/*.yaml`.
You will at least require the snakemake environment.

```
conda env create -f envs/snakemake.yaml
```

### Environments for the demo

The [demo](quickstart.md) (preprocessing, integration and metrics on the shipped example dataset) requires the following environments:

```
for env in snakemake scanpy bbknn scanorama scib plots funkyheatmap; do
    bash envs/install_environment.sh -f envs/$env.yaml
done
```

`envs/install_environment.sh` creates an environment if it does not exist yet and updates it otherwise.

The [funkyheatmap](https://funkyheatmap.github.io/funkyheatmap/) R package used for the summary plots of the `metrics` module is not available on conda.
Install it from CRAN into the `funkyheatmap` environment:

```
conda run -n funkyheatmap Rscript -e "install.packages('funkyheatmap', repos='https://cloud.r-project.org')"
```

Otherwise, the pipeline installs it on first use, which can fail when two plotting jobs try to install it at the same time (rerunning the pipeline resolves this).

**Typical install time:** ~10 minutes for the demo environments (measured on a laptop with Apple M1 Pro and an empty conda package cache, 0.5–1.5 minutes per environment).

> 💡 **Tip** On Apple Silicon, the `funkyheatmap` environment must be installed under emulation, because not all of its R dependencies are available for `osx-arm64`:
>
> ```
> CONDA_SUBDIR=osx-64 conda env create -f envs/funkyheatmap.yaml
> ```
>
> See [Troubleshooting](../principles/troubleshooting.md#working-on-apple-silicon).

### All environments

You can install the other environments as needed, for different parts of the workflow (modules, rules).
To install all environments at once, run:

```
bash envs/install_all_environments.sh
```

Use `-c mamba` to install with mamba and `-n` for a dry run that only lists what would be installed.
Note that GPU environments (`rapids_singlecell`, `scvi-tools`, `scarches`) require a machine with an NVIDIA GPU, see [Working with GPUs](../principles/troubleshooting.md#working-with-gpus).
