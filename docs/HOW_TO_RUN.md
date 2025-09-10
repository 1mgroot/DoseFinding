# How to Run the Simulation

There are two primary ways to run the Bayesian trial simulation:

1.  **Interactive Simulation with the Quarto Notebook (Recommended)**
2.  **Running the R Script Directly**

---

### 1. Interactive Simulation with the Quarto Notebook

This is the recommended method for running the simulation, as it provides an interactive environment to configure, run, and visualize the results.

**Steps:**

1.  **Open `notebooks/simulation_notebook.qmd`** in RStudio or your preferred R environment.
2.  **Modify the Configuration:** In the "Configuration" code chunk, you can adjust the `trial_config` list and the data simulation parameters (`p_YI`, `p_YT_given_I`, `p_YE_given_I`, etc.) to match your desired scenario.
3.  **Run the Simulation:** Run the code chunks in the notebook sequentially. The notebook will automatically source the necessary files from the `src/` directory, run the simulation, and generate plots and summary tables.

### 2. Running the R Script Directly

You can also run the simulation by executing the `main.R` script directly. This is useful for quick tests or for running the simulation in a non-interactive environment.

**Steps:**

1.  **Open `src/core/main.R`** in your R console or RStudio.
2.  **Source the Script:** Run the following command in your R console:

    ```R
    source("src/core/main.R")
    ```

This will execute the `run_trial_simulation` function with the default parameters defined in `src/core/config.R` and save the plots to the `results/plots` directory.
