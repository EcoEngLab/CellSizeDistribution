# CellSizeDistribution

This project provides tools for analyzing microbial cell size distributions by integrating cell morphology data from BacDive with amplicon sequencing data, and performing statistical and visualization analyses. The repository includes scripts for retrieving, calculating, merging, plotting, and statistically testing cell size data for bacteria and archaea.

## Features

- Retrieve cell morphology information (length, width, shape) from BacDive
- Calculate minimum, maximum, average, and geometric mean volumes and surface areas from cell size data
- Merge cell size database with amplicon sequencing data at the family or genus level
- Plot histograms and scatterplots of cell sizes by family or genus
- Perform permutation tests on cell size distributions
- Handle CSV import/export for all data processing steps

### Requirements

- Python 3.8+
- R (for R scripts)
- Git
- VS Code (recommended)
- VS Code Python extension
- VS Code R extension (for R scripts)
- VS Code PDF viewer extension (for viewing plots)

## Usage

- Run the main scripts in the `code/` directory according to your analysis needs:
  - Retrieve cell size data: `RetrieveCellSizeBacdive.py`
  - Calculate cell volumes: `VolumeCalc.Jolan.py`
  - Merge with amplicon data: `MergingOnFamilyLevel.py`, `MergingOnGenusLevel.py`
  - Plot histograms: `VolumeOfCellsByFamily.Bacteria.py`
  - Plot scatterplots: `FamilyAIC.R`
  - Run permutation tests: `ShufflingFamily.Bacteria`, `ShufflingGenus.Bacteria`

- You can add or remove more species to the relevant CSVs. Ensure header names stay the same.

## Contact/Support

For questions or support, please email ching.wan24@imperial.ac.uk
