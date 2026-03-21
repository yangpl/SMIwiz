# SMIwiz

SMIwiz is a C/MPI toolbox for acoustic seismic modelling, imaging, and inversion. It uses high-order staggered-grid finite-difference time-domain propagation and includes workflows for forward modelling, RTM, LSRTM, FWI, source estimation, migration deconvolution, and wavefield decomposition.

The repository contains the core solver, a small GTK launcher, reproducible synthetic examples, and reference materials for the published methods.

## Scope

SMIwiz currently supports:

- 2D and 3D acoustic modelling on staggered grids
- Isotropic, VTI, and TTI model input paths
- Forward modelling and data generation
- Full waveform inversion (FWI) and gradient evaluation
- Reverse time migration (RTM)
- Data-domain least-squares RTM (LSRTM)
- Source inversion
- Angle-domain common-image-gather (ADCIG) extraction
- Point-spread-function (PSF) Hessian calculation
- Migration deconvolution with PSF-based PCGNR or FFT-Wiener filtering
- Up-going and down-going wavefield separation

## Publications

If you use this code in research, cite the relevant paper(s):

- Pengliang Yang, “SMIwiz: An integrated toolbox for multidimensional seismic modelling and imaging”, Computer Physics Communications 295, 109011, 2024. DOI: [10.1016/j.cpc.2023.109011](https://doi.org/10.1016/j.cpc.2023.109011)
- Zhengyu Ji and Pengliang Yang, “SMIwiz-2.0: Extended functionalities for wavefield decomposition, linearized and nonlinear inversion”, Computer Physics Communications 309, 109503, 2025. DOI: [10.1016/j.cpc.2025.109503](https://doi.org/10.1016/j.cpc.2025.109503)
- Pengliang Yang and Zhengyu Ji, “A comparative study of data- and image-domain LSRTM under velocity-impedance parametrization”, Computers & Geosciences 208, 106091, 2026. DOI: [10.1016/j.cageo.2025.106091](https://doi.org/10.1016/j.cageo.2025.106091)

Project talk: <https://cassyni.com/events/P4W1QfiGXffuSf6Rv5VzJZ>

## Repository Layout

- `src/`: core C implementation and `Makefile`
- `include/`: public headers and shared data structures
- `bin/`: output directory for compiled executables
- `doc/`: papers and supporting technical references
- `GUI/`: GTK3 launcher for editing `inputpar.txt` and running jobs
- `run_fwd/`: minimal forward-modelling example
- `run_fwi2d/`: small 2D FWI example in a layered model
- `run_marmousi/`: 2D Marmousi FWI example
- `run_rtm/`: 2D RTM example
- `run_lsrtm2d/`: 2D two-parameter LSRTM example
- `run_lsrtm3d/`: 3D two-parameter LSRTM example
- `run_decon_1par/`: one-parameter migration deconvolution workflows
- `run_updown/`: up/down wavefield separation example
- `run_fbrec3d/`: 3D wavefield reconstruction example
- `run_fwi3d/`: 3D FWI example assets and survey generation
- `run_fwd_aniso/`: anisotropic forward-modelling example
- `run_Viking/`, `run_Viking_85shots/`: Viking Graben real-data workflows
- `run_rwi_grad/`: gradient-generation and related inversion experiments

## Requirements

- Linux
- MPI compiler and runtime, typically OpenMPI or MPICH
- FFTW3 development libraries
- `make`
- Optional: `gfortran` for helper programs in some example directories
- Optional: GTK3 development packages for the GUI

## Build

The main build uses `mpicc` and expects the FFTW installation prefix through the `fftw3` environment variable.

```bash
export fftw3=/path/to/fftw3
cd src
make
```

This creates the solver at `bin/SMIwiz`.

If your MPI wrapper or FFTW location differs, adjust [`src/Makefile`](/home/pyang/Documents/SMIwiz/src/Makefile).

## Running the Solver

SMIwiz is launched with MPI and a flat parameter list, usually generated from an `inputpar.txt` file:

```bash
mpirun -n <ranks> ../bin/SMIwiz $(cat inputpar.txt)
```

Most example directories already contain a `run.sh` or similarly named script that writes `inputpar.txt` and starts a job.

### Main `mode` Values

The executable entry point in [`src/main.c`](/home/pyang/Documents/SMIwiz/src/main.c) dispatches work by `mode=`:

- `0`: forward modelling
- `1`: time-domain FWI
- `2`: RTM
- `3`: data-domain LSRTM
- `4`: FWI gradient building
- `5`: source inversion
- `6`: ADCIG extraction
- `7`: PSF Hessian computation
- `8`: iterative migration deconvolution via PSF
- `9`: migration deconvolution via FFT-Wiener filtering
- `10`: up/down wavefield separation

### Common Required Inputs

Most runs need at least:

- grid and time parameters: `n1`, `n2`, optional `n3`, `d1`, `d2`, optional `d3`, `nt`, `dt`
- acquisition description: `acquifile=...`
- model files: `vpfile=...`, `rhofile=...`
- source wavelet: `stffile=...` unless `mode=5`
- inversion or imaging options depending on the selected mode

Model files are read as binary `float` arrays with dimensions matching the declared grid size.

## Quick Start

### 1. Build the executable

```bash
export fftw3=/path/to/fftw3
cd src
make
```

### 2. Run the small 2D FWI example

```bash
cd ../run_fwi2d
bash run.sh
```

That script writes `inputpar.txt` and launches:

```bash
mpirun -n 2 ../bin/SMIwiz $(cat inputpar.txt)
```

### 3. Run the 2D RTM example

```bash
cd ../run_rtm
bash run_rtm.sh
```

The RTM example is configured for a larger MPI job in the provided script, so adjust the `mpirun -n` value if your machine has fewer ranks available.

## Example Notes

- `run_fwi2d/` is the easiest place to start. It is small, self-contained, and shows the standard `inputpar.txt` workflow clearly.
- `run_marmousi/` and `run_fwd_aniso/` are useful for more realistic 2D tests and anisotropic inputs.
- `run_decon_1par/` expects a staged workflow: first compute the PSF (`mode=7`), then run deconvolution with `mode=8` or `mode=9` depending on method.
- `run_fbrec3d/` requires generating `acqui.txt` locally before running because the full file is too large to keep in the repository.
- `run_Viking/` is the real-data workflow. Follow [`run_Viking/0_readme.txt`](/home/pyang/Documents/SMIwiz/run_Viking/0_readme.txt) for its step-by-step process.

## GUI

The GTK launcher in `GUI/` provides a simple interface for editing `inputpar.txt`, selecting the run directory, and launching `mpirun`.

Build and run it with:

```bash
cd GUI
make
./gui
```

See [`GUI/README.md`](/home/pyang/Documents/SMIwiz/GUI/README.md) for the current GUI behavior and requirements.

## Reproducibility and Version Notes

A few examples refer to specific historical revisions in their local notes. If you are trying to reproduce figures from a paper exactly, check the run directory for version-specific instructions before changing scripts or parameters.

## License

This project is distributed under the GPL-3.0 license. See [`LICENSE`](/home/pyang/Documents/SMIwiz/LICENSE).
