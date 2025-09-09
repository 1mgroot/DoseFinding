## How to Run the Simulation

This document explains how to run the Bayesian trial simulation.

### 1. Prerequisites

Make sure you have R installed on your system. You will also need the following R libraries:

- `dplyr`
- `tidyr`
- `isotone`
- `purrr`
- `ggplot2`
- `Iso`
- `testthat`

You can install them using the following command in your R console:

```R
install.packages(c("dplyr", "tidyr", "isotone", "purrr", "ggplot2", "Iso", "testthat"))
```

### 2. Running the Simulation

To run the simulation, open your R console and execute the following command:

```R
source("main.R")
```

This will run the entire simulation pipeline, from data simulation to model fitting and visualization. The results, including the selected optimal dose (OD) and plots, will be saved in the `results` directory.

### 3. Configuration

The simulation parameters can be modified in the `config.R` file. This includes the dose levels, number of stages, cohort size, and decision thresholds.

### 4. Testing

To run the unit tests, execute the following command in your R console:

```R
source("test_main.R")
```

This will run the tests for the main module and print the results to the console.
