# OpenPlumes

OpenPlumes is an open-source research project for three-dimensional groundwater contaminant plume modeling and visualization.

The long-term goal is a standalone scientific engine that can support Python scripts, Jupyter notebooks, QGIS, web applications, Blender, and other geospatial interfaces. The current release is an experimental Python prototype focused on a transparent, reproducible modeling workflow.

## Current prototype

The standalone workflow currently supports:

- loading synthetic monitoring-well data from CSV;
- projecting longitude and latitude into a metric coordinate reference system;
- creating three-dimensional well points;
- reading a polygon model boundary;
- generating a three-dimensional prediction grid;
- fitting a SciPy radial basis function interpolator;
- predicting contaminant concentrations throughout the grid;
- extracting concentration isosurfaces;
- generating nearest-neighbor concentration profiles;
- and displaying the results interactively with PyVista.

## Installation

Open a terminal in the repository directory and create the Conda environment:

```bash
conda env create -f environment.yaml
conda activate openplumes
```

If the `openplumes` environment already exists, activate it directly:

```bash
conda activate openplumes
```

## Run the example

From the repository root, run:

```bash
python main.py
```

The example uses the synthetic well data and model boundary stored in `Examples/`. A PyVista window should open containing the wells and interpolated plume isosurfaces.

The example data are entirely synthetic and do not represent observed contamination at a real location.

## Current project layout

```text
OpenPlumes/
├── engine/             Standalone scientific functions
├── Examples/           Synthetic example inputs
├── Notebooks/          Research and algorithm-development notebooks
├── openplumes/         Earlier experimental QGIS plugin
├── environment.yaml    Conda environment definition
└── main.py             Standalone example workflow
```

The `engine/` directory and `main.py` contain the current standalone prototype. The lowercase `openplumes/` directory contains the earlier QGIS plugin experiment. QGIS is not required to run `main.py`. The plugin will eventually become a lightweight interface to the standalone engine.

## Scientific status

OpenPlumes is currently a research preview, not a validated environmental decision-making tool.

Known limitations include:

- RBF interpolation can overshoot the observed concentration range;
- predicted concentrations may become negative;
- predictions outside the monitoring-well network are extrapolations;
- disconnected or boundary-clipped isosurfaces may be interpolation artifacts;
- hydrogeologic anisotropy is not yet represented explicitly;
- and sample depth currently assumes a flat reference surface unless elevation data are supplied.

These limitations are active research topics. Scientific correctness, diagnostics, and validation against trusted libraries take priority over optimization and interface development.

## Development approach

Jupyter notebooks are used as the research laboratory for understanding and validating algorithms. Functions are moved into `engine/` once their inputs, outputs, assumptions, and failure modes are sufficiently understood. The standalone `main.py` workflow verifies that the pipeline can run without hidden notebook state.

## Interactive examples

- [OpenPlumes interactive plume example 1](https://jmm688.github.io/scene-export%20(11).html)
- [OpenPlumes interactive plume example 2](https://jmm688.github.io/scene-export%20(17).html)

## License

See [LICENSE](LICENSE).
