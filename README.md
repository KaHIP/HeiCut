HeiCut
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.16+-064F8C.svg)](https://cmake.org/)
[![Linux](https://img.shields.io/badge/Linux-supported-success.svg)](https://github.com/KaHIP/HeiCut)
[![GitHub Stars](https://img.shields.io/github/stars/KaHIP/HeiCut)](https://github.com/KaHIP/HeiCut/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/KaHIP/HeiCut)](https://github.com/KaHIP/HeiCut/issues)
[![Last Commit](https://img.shields.io/github/last-commit/KaHIP/HeiCut)](https://github.com/KaHIP/HeiCut/commits)
[![ALENEX'26](https://img.shields.io/badge/ALENEX'26-10.1137/1.9781611978957.13-blue)](https://doi.org/10.1137/1.9781611978957.13)
[![arXiv](https://img.shields.io/badge/arXiv-2504.19842-b31b1b.svg)](https://arxiv.org/abs/2504.19842)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17140472.svg)](https://doi.org/10.5281/zenodo.17140472)
[![Heidelberg University](https://img.shields.io/badge/Heidelberg-University-c1002a)](https://www.uni-heidelberg.de)
=====

<p align="center">
  <img src="https://raw.githubusercontent.com/KaHIP/HeiCut/main/logo/banner.png" alt="HeiCut Logo" width="900"/>
</p>

**HeiCut** is a highly efficient, exact solver for the **minimum cut problem in hypergraphs** using FPT kernelization. Given a hypergraph, the minimum cut problem asks for a partition of the vertices into two non-empty sets such that the total weight of hyperedges crossing the partition is minimized. HeiCut performs repeated rounds of provably exact reduction rules that preserve the minimum cut while drastically shrinking instance size, then applies an exact solver to the reduced hypergraph. Part of the [KaHIP](https://github.com/KaHIP) organization.

| | |
|:--|:--|
| **What it solves** | Exact minimum cut in weighted and unweighted hypergraphs |
| **Techniques** | FPT kernelization, provably exact reduction rules, label propagation coarsening |
| **Solvers** | HeiCut (kernelization + exact solve), Relaxed BIP (Gurobi ILP), Trimmer, vertex-ordering solvers |
| **Requires** | C++17 compiler (GCC 7+), CMake 3.16+, Boost, oneTBB, hwloc, SparseHash, Gurobi, Mt-KaHyPar |

## Quick Start

### Build from source

```bash
git clone https://github.com/KaHIP/HeiCut.git && cd HeiCut

# Install dependencies (Ubuntu)
sudo apt-get install libtbb-dev libhwloc-dev libboost-program-options-dev

# Install Mt-KaHyPar
./install_mtkahypar.sh

# Build
mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

### Run

```bash
# Run HeiCut with tight ordering
./build/kernelizer path/to/hypergraph.hgr --ordering_type=tight

# With label propagation
./build/kernelizer path/to/hypergraph.hgr --ordering_type=tight --lp_num_iterations=1

# ILP solver (2-hour timeout)
./build/ilp path/to/hypergraph.hgr --ilp_timeout=7200
```

---

## How It Works

The minimum cut problem in hypergraphs is NP-hard. HeiCut applies **fixed-parameter tractable (FPT) kernelization** to reduce the hypergraph before solving. The data reduction rules provably preserve the exact minimum cut value while shrinking the hypergraph, often dramatically. The reduced hypergraph is then solved using an exact solver.

The pipeline:
1. **Load** a hypergraph (hMETIS or METIS format)
2. **Reduce** via provably exact kernelization rules, optionally with label propagation coarsening
3. **Solve** the reduced hypergraph exactly (ILP, vertex-ordering, or trimmer)
4. **Extract** the minimum cut value for the original hypergraph

---

## Executables

| Binary | Description |
|:-------|:------------|
| `kernelizer` | Main HeiCut solver (sequential) |
| `kernelizer_parallel` | Parallel HeiCut solver |
| `ilp` | Relaxed BIP solver (Gurobi) |
| `ilp_parallel` | Parallel ILP solver |
| `trimmer` | Trimmer algorithm |
| `submodular` | Vertex-ordering solver |
| `submodular_parallel` | Parallel vertex-ordering solver |
| `dumbbell_generator` | Synthetic dumbbell hypergraph generator |
| `kcore_generator` | (k,2)-core benchmark generator |

All executables support `--help` to list available arguments.

---

## Command Line Usage

### HeiCut (kernelizer)

```bash
# Tight ordering, no label propagation
./build/kernelizer PATH_TO_HYPERGRAPH --ordering_type=tight

# Tight ordering + 1 round of label propagation
./build/kernelizer PATH_TO_HYPERGRAPH --ordering_type=tight --lp_num_iterations=1

# Parallel execution
./build/kernelizer_parallel PATH_TO_HYPERGRAPH --ordering_type=tight
```

> **Tip:** Use `--verbose` to view detailed reduction performance.

### Relaxed BIP (ILP)

```bash
# ILP solver with 2-hour timeout
./build/ilp PATH_TO_HYPERGRAPH --ilp_timeout=7200

# Parallel ILP solver
./build/ilp_parallel PATH_TO_HYPERGRAPH --ilp_timeout=7200
```

### Trimmer

```bash
./build/trimmer PATH_TO_HYPERGRAPH --ordering_type=tight
```

### Vertex-Ordering Solver

```bash
./build/submodular PATH_TO_HYPERGRAPH --ordering_type=tight

# Parallel
./build/submodular_parallel PATH_TO_HYPERGRAPH --ordering_type=tight
```

### Generators

```bash
# Generate a (k,2)-core hypergraph
./build/kcore_generator PATH_TO_HYPERGRAPH PATH_TO_OUTPUT

# Generate a dumbbell hypergraph
./build/dumbbell_generator PATH_TO_OUTPUT
```

---

## Hypergraph Formats

Two input formats are supported. The default is **hMETIS** (see [hMETIS manual](https://course.ece.cmu.edu/~ee760/760docs/hMetisManual.pdf)). **METIS** format is also supported.

**hMETIS format** (unweighted):
```
% comment line (optional)
num_hyperedges num_vertices
vertex_1 vertex_2 ...
vertex_3 vertex_4 vertex_5 ...
...
```

**hMETIS format** (weighted, with `1` flag):
```
% comment line (optional)
num_hyperedges num_vertices 1
weight_1 vertex_1 vertex_2 ...
weight_2 vertex_3 vertex_4 vertex_5 ...
...
```

Example files are included in the `examples/` directory.

---

## Building from Source

### Dependencies

- A 64-bit **Linux** operating system
- A modern **C++17** compiler (`g++` >= 7 recommended)
- [CMake](http://www.cmake.org/) (>= 3.16)
- [Boost.Program_options](http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html) (>= 1.69)
- [oneTBB](https://github.com/oneapi-src/oneTBB) (>= 2021.5.0)
- [hwloc](https://www.open-mpi.org/projects/hwloc/)
- [SparseHash](https://github.com/sparsehash/sparsehash)
- [Gurobi](https://www.gurobi.com/) (used as LP solver)
- [Mt-KaHyPar](https://github.com/kahypar/mt-kahypar/tree/0ef674ad44c35fb4f601a7eddd3f4f23f0d5d60a), commit `0ef674a`

### Install core dependencies (Ubuntu/Debian)

```bash
sudo apt-get install libtbb-dev libhwloc-dev libboost-program-options-dev
```

> **Note:** `libtbb-dev` may be outdated on some distributions. If so, clone and build [oneTBB](https://github.com/oneapi-src/oneTBB) locally:

```bash
git clone https://github.com/oneapi-src/oneTBB.git
cd oneTBB && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr -DTBB_TEST=OFF ..
sudo cmake --build .
sudo cmake --install .
```

### Install SparseHash

```bash
git clone https://github.com/sparsehash/sparsehash
cd sparsehash
./configure
make install
```

### Install Gurobi

HeiCut uses **Gurobi** as its LP solver.

1. Create a [Gurobi account](https://www.gurobi.com/).
2. Download [Gurobi for Linux](https://www.gurobi.com/downloads/gurobi-software/).
3. Follow the [installation guide](https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer).
4. Obtain a license. Free academic licenses are available via the [official guide](https://www.gurobi.com/academia/academic-program-and-licenses/).
5. Place your `gurobi.lic` license file in the installation folder (e.g., `/opt/gurobi1203/`).
6. Add environment variables (adjust the path to match your version):

```bash
export GUROBI_HOME="/opt/gurobi1203/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```

### Install Mt-KaHyPar

HeiCut depends on **Mt-KaHyPar** (commit `0ef674a`).

To install automatically:

```bash
./install_mtkahypar.sh
```

This builds the library and places `libmtkahypar.so` in `HeiCut/extern/mt-kahypar-library/`.

Manual build instructions (if preferred):

```bash
git clone --depth=2 --recursive https://github.com/kahypar/mt-kahypar.git
cd mt-kahypar
git checkout 0ef674a
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make install.mtkahypar   # may need sudo
```

**Note:** If installed locally, the build will exit with an error due to missing permissions. However, the library is still built successfully and is available in the build folder. Locate `libmtkahypar.so` (usually in `build/lib/`) and copy it into `HeiCut/extern/mt-kahypar-library/`.

#### Known issue (rare): `growt` ref/_mref compilation error

A few users have encountered a compilation error in Mt-KaHyPar's `growt` dependency:

```
.../external_tools/growt/data-structures/migration_table_iterator.hpp:68:22: error:
'... migration_table_mapped_reference ...' has no member named 'ref'; did you mean '_mref'?
68 |                 sref.ref.refresh();
|                      ^~~
|                      _mref
```

If you see this and the program does not compile, apply the following manual fix in the `growt` source and rebuild:

1. Open:
```
<your-mtkahypar-source>/external_tools/growt/data-structures/migration_table_iterator.hpp
```

2. In the class `migration_table_mapped_reference`, edit the constructor's initializer list. Change (possibly line 57):
```cpp
: _tab(table), _version(ver), _mref(mref)
```
to:
```cpp
: _tab(table), _version(ver), _mref(mref), ref(_mref)
```

3. Insert (possibly line 133):
```cpp
public:
   base_mapped_reference& ref;
```

4. Rebuild Mt-KaHyPar (or rerun your previous build command).

### Build HeiCut

After installing all dependencies:

```bash
git clone https://github.com/KaHIP/HeiCut.git && cd HeiCut
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Binaries are placed in `HeiCut/build/`.

---

## Reproducing Paper Results

Scripts to reproduce all experimental results from the paper are located in `experiments/`, organized per dataset:

- `medium_weighted`, `medium_unweighted`
- `large_weighted`, `large_unweighted`
- `k-core_weighted`, `k-core_unweighted`

### Benchmark Datasets

We evaluate on three datasets:

- **M<sub>HG</sub>** (488 medium instances): weighted and unweighted
- **L<sub>HG</sub>** (94 large instances, up to 139M vertices): weighted and unweighted
- **(k,2)-core** (44 synthetic instances with non-trivial cuts): weighted and unweighted

Download all datasets here:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17142170.svg)](https://doi.org/10.5281/zenodo.17142170)

Extract into the repository root (`/HeiCut/`):

- `/HeiCut/med_set/`
- `/HeiCut/large_set/`
- `/HeiCut/k,2-core_benchmark/`

**Note:** To save time and disk space, only weighted versions of the hypergraphs are provided. For reproducing results on unweighted versions, pass `--unweighted` to each algorithm to process the weighted hypergraphs as unweighted (all edge weights set to 1).

### Running Experiments

Each dataset folder contains a master script (e.g., `perform_medium_weighted_experiments.sh`). A global script `perform_all_experiments.sh` is also provided, though running all experiments at once is very time-consuming.

**Important:** Avoid using special characters (spaces, #, %, &) in the HeiCut directory path, as they can break the experimental and plotting scripts.

**Experiment dependencies:**

```bash
sudo apt install parallel time
```

**Setup notes:**
- Time limits: 2 hours per algorithm
- Memory limits: 100 GB for M<sub>HG</sub>, 300 GB for L<sub>HG</sub> and (k,2)-core
- By default, scripts run 1 instance at a time (adjustable for machines with more memory)

### Results Output

Results are written to `experiments/<dataset>/generated/all_results/`. Each algorithm produces CSV summaries and per-instance results. In the output `all_results.csv`, each row contains statistics for an instance. The first three columns correspond to minimum cut, time, and memory. If an algorithm fails on an instance, the minimum cut column is blank. All algorithms except Relaxed BIP return the exact minimum cut if successful.

### Plotting

Each dataset folder contains a `plot/` subfolder. Once experiments are complete:

```bash
cd experiments/medium_weighted/plot
./plot_all.sh
```

This generates performance profile plots comparing all algorithms on memory usage, runtime, and minimum cut, replicating the paper's figures.

**Plotting dependencies:**
- **R** (tested with R 4.3+)
- **LaTeX** with TikZ support (a full texlive installation is recommended)

Required R packages:
```
ggplot2, plyr, dplyr, RColorBrewer, tikzDevice,
gridExtra, egg, ggpubr, stringr, stringi, ggrepel
```

The R scripts will automatically attempt to install missing packages. If installation fails, install them manually:
```r
install.packages(c(
  "ggplot2","plyr","dplyr","RColorBrewer","tikzDevice",
  "gridExtra","egg","ggpubr","stringr","stringi","ggrepel"
), repos = "https://cloud.r-project.org")
```

### Custom Experiments

Configurable scripts are available in `experiments/`:

- `generate_experiments.sh` -- select algorithms and parameters
- `run_experiments.sh` -- set time limits and parallelism
- `extract_results.sh` -- collect results

Provide hypergraph paths in `experiments/hypergraphs.txt`.

---

## Related Projects

| Project | Description |
|:--------|:------------|
| [KaHIP](https://github.com/KaHIP/KaHIP) | Karlsruhe High Quality Graph Partitioning |
| [VieCut](https://github.com/KaHIP/VieCut) | Shared-memory parallel minimum cut algorithms |
| [fpt-max-cut](https://github.com/KaHIP/fpt-max-cut) | FPT kernelization for the maximum cut problem |
| [Mt-KaHyPar](https://github.com/kahypar/mt-kahypar) | Multi-threaded Karlsruhe Hypergraph Partitioner |

---

## Licence

HeiCut is free software provided under the MIT License.
If you publish results using our algorithms, please cite:

```bibtex
@inproceedings{DBLP:conf/alenex/Chhabra0UW26,
  author       = {Adil Chhabra and
                  Christian Schulz and
                  Bora U{\c{c}}ar and
                  Loris Wilwert},
  editor       = {Rezaul Chowdhury and
                  Simon J. Puglisi and
                  Bin Ren and
                  Nate Veldt},
  title        = {Exact Minimum Cuts in Hypergraphs at Scale},
  booktitle    = {Proceedings of the 28th Symposium on Algorithm Engineering and Experiments,
                  {ALENEX} 2026, Vancouver, BC, Canada, January 11-12, 2026},
  pages        = {169--181},
  publisher    = {{SIAM}},
  year         = {2026},
  url          = {https://doi.org/10.1137/1.9781611978957.13},
  doi          = {10.1137/1.9781611978957.13},
}
```

The source code and benchmark datasets are permanently archived on [Zenodo (Software)](https://doi.org/10.5281/zenodo.17140472) and [Zenodo (Dataset)](https://doi.org/10.5281/zenodo.17142170).

## References

1. C. J. Alpert, *The ISPD98 Circuit Benchmark Suite*, ISPD 1998. [DOI](https://doi.org/10.1145/274535.274546)
2. N. Viswanathan et al., *The DAC 2012 Routability-Driven Placement Contest*, DAC 2012. [DOI](https://doi.org/10.1145/2228360.2228500)
3. T. A. Davis and Y. Hu, *The SuiteSparse Matrix Collection*, ACM TOMS 2011. [DOI](https://doi.org/10.1145/2049662.2049663)
4. A. Belov et al., *The SAT Competition 2014*. [Link](https://satisfiability.org/competition/2014/)
