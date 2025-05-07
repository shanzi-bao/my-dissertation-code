# my-dissertation-code

This repository contains the R code associated with my **double-semester project dissertation**.  
Each script or folder corresponds to a specific section of the dissertation and can be run directly after downloading the repository.

> ✅ All code has been developed and tested under **R version 4.4.2**.  
> For best compatibility, it is recommended to run the scripts using this version.

---

## 📁 File Overview

- **GEV distribution plot.R**  
  → Generates plots of the GEV distribution and its tail behavior  
  → Corresponds to Section **2.1.1** of the dissertation

- **Smith model spatial dependence structure picture.R**  
  → Visualizes the dependence structure of the Smith max-stable model  
  → Corresponds to Section **2.1.4**

- **Different ABC algorithm comparisons.R**  
  → Compares multiple ABC algorithms in a simulation study  
  → Corresponds to Section **2.3.5**

- **ABC regression adjustment comparison.R**  
  → Demonstrates the effect of regression adjustment in ABC postprocessing  
  → Corresponds to Section **2.3.6**

- **grid example to present different batch strategy.R**  
  → Visualization of different batch strategies in a grid  
  → Corresponds to Section **3.2.3**

---

## 🧪 Experiments

### 📁 `Experiment1/` – Numerical Experiment I (Section 4.4)

This experiment focuses on a **regular spatial grid** setup.  
To reproduce the results:

- Run `rnum.R` for **MCS(NUM)-ABC**
- Run `rcs.R` for **ICS-ABC**

Both use the **adaptive ABC-PMC** algorithm.  
The default setup assumes **k = 2**; if you'd like to try other values, simply replace `2` with your desired `k` inside the scripts.

Supporting files:

- `observed_data.R`: generates the observed dataset  
- `supplementary_function.R`: utility functions required by both models  
- `summary_statistics_calculation.cpp`: C++ code used for summary statistic computation

---

### 📁 `Experiment2/` – Numerical Experiment II (Section 4.5)

This experiment is based on a **realistic station network**.  
To reproduce the results:

- Run `num.R` for **MCS(NUM)-ABC**
- Run `cs.R` for **ICS-ABC**

All instructions are the same as in Experiment 1.  
Again, **k = 2** is the default; modify to change the dimension as needed.

---


## 🧠 ABC Regression Adjustment (Custom for Section 4.4 & 4.5)

Although the `abc` package includes a basic regression adjustment option,  
our experimental models (especially in Sections **4.4** and **4.5**) involve complex structured summary statistics,  
which can not be directly supported by standard postprocessing tools.

Therefore, we implemented a custom function `abc_regression_adjustment_v2()` that uses the local linear regression.

This function is used **only** to postprocess results in the two numerical experiments and is not needed elsewhere in the project.



## ✅ How to Use

1. Clone or download the entire repository.
2. Place all files under the same R project directory.
3. Open the `.Rproj` file in **RStudio**.
4. Run any script directly (e.g. `rnum.R`, `num.R`) to reproduce the corresponding result.
5. Ensure your R environment is version **4.4.2** for full compatibility.

---

## 📌 Note

- All scripts are self-contained and reproducible.
- Figures generated match the ones shown in the dissertation.
- Adaptive ABC-PMC is used in all experiments unless otherwise stated.

---

## 📄 License

This repository is for academic purposes only and is shared to accompany my **double-semester project**.

## 📚 References

