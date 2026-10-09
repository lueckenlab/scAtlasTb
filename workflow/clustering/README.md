# Clustering

This module provides a workflow to perform clustering at multiple resolutions and at multiple subclustering levels (iterative sub-clustering).
The code is optimised for memory efficiency and speed, allowing for clustering of large datasets.

## Configuration

Example configuration for clustering can be found in `workflow/clustering/test/config.yaml` showcasing the setup for hierarchical (iterative) sub-clustering.

```yaml
DATASETS:

  test_level1_multiresolution:
    input:
      clustering: test/input/pbmc68k.h5ad
    clustering:
      recompute_umap: true
      recompute_neighbors: true
      neighbors:
        n_neighbors: 30
        use_rep: X_pca
      algorithm:
        - leiden
        # - louvain
      resolutions:
        - 0.2
        - 0.3
        - 0.4
      umap_colors:
        - bulk_labels
        - batch
        
  test_hierarchical:
    input:
      clustering: test/input/pbmc68k.h5ad
    clustering:
      neighbors_key: neighbors
      algorithm: leiden
      resolutions:
        - 0.5
      hierarchy:
        1: 0.1
        3: 0.2
      umap_colors:
        - bulk_labels
        - batch
        - n_genes
        
  test_hierarchical_int:
    input:
      clustering: test/input/pbmc68k.h5ad
    clustering:
      neighbors_key: neighbors
      algorithm: leiden
      use_gpu: true
      resolutions:
        - 0.5
      hierarchy: 3
      umap_colors:
        - bulk_labels
        - batch
```

## Parameters

### Clustering Parameters

- **`algorithm`**: Clustering algorithm to use
  - Single algorithm: `leiden` or `louvain`
  - Multiple algorithms: List of algorithms (e.g., `[leiden, louvain]`)

- **`resolutions`**: List of clustering resolutions to test
  - Example: `[0.2, 0.4, 0.8, 1.0, 1.6]`

- **`neighbors_key`**: Key for pre-computed neighbors in `adata.uns`
  - Default: Uses existing neighbors if available

- **`hierarchy`**: Configuration for iterative sub-clustering
  - Dictionary format: `{level: resolution}` (e.g., `{1: 0.1, 3: 0.2}`); levels up to the largest key are computed, levels missing from the dictionary use the value from `resolutions`
  - Integer format: Maximum hierarchy depth (e.g., `3`)

### Neighbor Computation

- **`recompute_neighbors`**: Boolean flag to recompute neighbors
  - `true`: Force recomputation of neighbors
  - `false` (default): Use existing neighbors if available

- **`neighbors`**: Parameters for neighbor computation when recomputing
  - `n_neighbors`: Number of neighbors (e.g., `30`)
  - `use_rep`: Representation to use (e.g., `X_pca`)

### UMAP Configuration

- **`recompute_umap`**: Boolean flag to recompute UMAP embedding
  - `true`: Force recomputation of UMAP
  - `false` (default): Use existing UMAP if available

- **`umap_colors`**: List of metadata columns to use for UMAP coloring
  - Example: `[bulk_labels, batch, n_genes]`

### Performance Options

- **`use_gpu`**: Boolean flag to enable GPU acceleration
  - `true`: Use GPU for clustering when available
  - `false` (default): Use CPU only

- **`n_cell_cpu`**: Threshold for forcing CPU computation - CPU is used when dataset has fewer than this many cells (default: 100000)

## Calling the clustering script from the command line

The clustering script is located at `workflow/clustering/scripts/clustering.py` and performs single-level clustering operations. For hierarchical clustering, you need to run the script multiple times in sequence.

The script automatically detects GPU availability and uses RAPIDS when possible for accelerated clustering. It falls back to scanpy implementation when needed.

Make sure you are using the `envs/scanpy.yaml` environment, which has all the necessary dependencies installed.

Check the help message for the script to see all available options:

```bash
python workflow/clustering/scripts/clustering.py -h
```

```
usage: clustering.py [-h] [--resolution RESOLUTION] [--algorithm {louvain,leiden}] [--level LEVEL] [--threads THREADS] [--overwrite] [--max_cluster_factor MAX_CLUSTER_FACTOR] [--clustering_args CLUSTERING_ARGS] [--neighbors_key NEIGHBORS_KEY] [--neighbors_args NEIGHBORS_ARGS] [--n_cell_cpu N_CELL_CPU]
                     input_file output_file

Run clustering on anndata file.

positional arguments:
  input_file            Input anndata file (.h5ad or .zarr)
  output_file           Output anndata file (.zarr)

options:
  -h, --help            show this help message and exit
  --resolution RESOLUTION
                        Clustering resolution (default: 1.0)
  --algorithm {louvain,leiden}
                        Clustering algorithm (default: leiden)
  --level LEVEL         Hierarchical clustering level (default: 1)
  --threads THREADS     Number of threads for parallel processing (default: 1)
  --overwrite           Overwrite existing clustering results
  --max_cluster_factor MAX_CLUSTER_FACTOR
                        Maximum cluster factor for heuristic check (default: 50)
  --clustering_args CLUSTERING_ARGS
                        Additional clustering arguments as JSON string
  --neighbors_key NEIGHBORS_KEY
                        Key for neighbors in adata.uns (default: neighbors)
  --neighbors_args NEIGHBORS_ARGS
                        Arguments for neighbors computation as JSON string
  --n_cell_cpu N_CELL_CPU
                        Force CPU when dataset has fewer than this many cells (default: 100000)
```

## Script Behavior

### Level 1 Clustering
For `--level 1`, the script:
1. Reads the input file with existing neighbors graph
2. Applies clustering using the specified algorithm and resolution
3. Outputs a zarr file containing only the clustering results in `.obs`

### Hierarchical Clustering (Level > 1)
For `--level > 1`, the script:
1. Reads clustering results from the previous level
2. Recomputes neighbors for each cluster separately using parallel processing
3. Applies sub-clustering within each cluster
4. Combines results with hierarchical naming (e.g., `cluster_subcluster`)

### GPU/CPU Handling
- Automatically detects NVIDIA GPU availability
- Uses RAPIDS implementation when GPU is available and beneficial
- Falls back to scanpy for small datasets (< `n_cell_cpu` cells)
- Includes a heuristic check on the number and size of GPU clusters, falling back to the scanpy (CPU) implementation if it fails

### Output Format
- Creates zarr-linked output preserving the original data structure
- Only writes the new clustering column to minimize storage
- Clustering results are named as `{algorithm}_{resolution}_{level}`

## Examples

### Basic clustering at level 1:
```bash
python workflow/clustering/scripts/clustering.py \
    input.h5ad \
    output_level1.zarr \
    --resolution 1.0 \
    --algorithm leiden
```

### Hierarchical clustering at level 2:
```bash
# First run level 1
python workflow/clustering/scripts/clustering.py \
    input.h5ad \
    output_level1.zarr \
    --resolution 1.0 \
    --level 1

# Then run level 2 (requires level 1 results)
python workflow/clustering/scripts/clustering.py \
    output_level1.zarr \
    output_level2.zarr \
    --resolution 0.5 \
    --level 2 \
    --threads 4 \
    --neighbors_args '{"n_neighbors": 15, "use_rep": "X_pca"}'
```

### Using custom clustering parameters:
```bash
python workflow/clustering/scripts/clustering.py \
    input.h5ad \
    output.zarr \
    --resolution 0.8 \
    --clustering_args '{"n_iterations": 5, "flavor": "igraph"}' \
    --overwrite
```

## Output Structure

Each run produces a single clustering result column:
- **Level 1**: `leiden_1.0_1` (algorithm_resolution_level)
- **Level 2**: `leiden_1.0_2` (inherits resolution from level 1)

To collect multiple clustering results, read and combine the `.obs` columns:

```python
import pandas as pd
from utils.io import read_anndata

# Combine results from multiple runs
obs_data = []
for file in ['output_level1.zarr', 'output_level2.zarr']:
    obs = read_anndata(file, obs='obs', verbose=False).obs
    obs_data.append(obs)

combined_obs = pd.concat(obs_data, axis=1)
```