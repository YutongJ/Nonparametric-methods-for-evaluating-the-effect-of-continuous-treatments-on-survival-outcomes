# README

This repository contains the code and supporting materials for the paper:  
**"A Class of Nonparametric Methods for Evaluating the Effect of Continuous Treatments on Survival Outcomes"**.

## Overview

In this work, we propose and implement a novel nonparametric hypothesis testing methods to evaluate the effect of continuous treatments (or exposures) on survival outcomes. This method address challenges presented with censored data, continuous exposure, and potential nonlinearity in the treatment-outcome relationship. 


The repository provides:

- Code for simulating datasets with continuous treatments and survival outcomes in three scenarios (no effect, monotone effect, and concave effect).
- Implementation of the proposed hypothesis testing.
- Reproducible scripts for the analyses and results presented in the paper.



## Features

1. **Simulation Framework**:
   - Generate survival data with a continuous treatment and two confounding factors.
   - Include censoring mechanisms and treatment-outcome relationships (none, monotone or concave).

2. **Hypothesis Testing Methods**:
   - Nonparametric hypothesis testing procedure described in the paper.
   - Nuisance parameter estimation:
   		- Density of the exposure: kernel smoothing method.
   		- Conditional probability of survival/censoring: survival super learner (ensemble learning algorithm).


## Repository Structure

```plaintext
├── Data_Application/
│   ├── R_func/								# Core helper functions
│   ├── step0_Application_preprocessing.R	# Preprocessing real-world datasets
│   ├── step1_Application_analysis.R		# Application of methods
│   └── step2_Application_visualization.R	# Results and visualization 
├── Simulations/
│   ├── R_func/								# Core helper functions
│   ├── step0_DGP_and_PSfit.R				# Data generation scripts and estimation of exposure density
│   ├── step1_simulation_clustering.R		# Reproducible scripts for simulations
│   ├── step2_simulation_outputs.R			# Organizing simulation outputs
│   ├── step3_simulation_summary.R			# Results summary
│   └── step4_simulation_visualization.R	# Results visualization
├── README.md                     # This file
└── LICENSE                       # License details
```

## Prerequisites

- **R (version >= 4.4.0)** with the following packages:
  - `survSuperLearner`
  - `CVXR`
  - `mgcv`
  - `hal9001`
  - `ggplot2`
  - `gridExtra`
  - `patchwork`
- Additional dependencies can be installed using the provided script: `setup/install_packages.R`.



## Getting Started

1. **Clone the Repository**:
	 ```bash
	git clone https://github.com/YutongJ/Nonparametric-methods-for-evaluating-the-effect-of-continuous-treatments-on-survival-outcomes.git
	cd Nonparametric-methods-for-evaluating-the-effect-of-continuous-treatments-on-survival-outcomes
	```

2. **Install Required Packages**:
	Run the package installation script:
	```R
	# List of required packages
	required_packages <- c(
	  "survSuperLearner",	# For survival analysis
	  "CVXR",				# For hypothesis testing procedure
	  "hal9001",			# For estimating density of exposure
	  "ggplot2",			# For data visualization
	  "gridExtra",			# For data visualization
	  "patchwork"			# For data visualization
	)
	# Function to check and install missing packages
	install_if_missing <- function(packages) {
	  installed_packages <- rownames(installed.packages())
	  for (pkg in packages) {
	    if (!pkg %in% installed_packages) {
	      message(paste("Installing package:", pkg))
	      install.packages(pkg, dependencies = TRUE)
	    } else {
	      message(paste("Package already installed:", pkg))
	    }
	  }
	}
	# Install required packages
	install_if_missing(required_packages)
	```

<!-- 3. **Run Example Analysis**:
	Execute the example script to demonstrate the methods:
	```R
	source(paste0(dir_fun, "DGP_complex_ceiling.R"))
	source(paste0(dir_fun, "EstPropScore.R"))
	source(paste0(dir_fun, "getIF.R"))
	source(paste0(dir_fun, "causalSS_SL.R"))
	source(paste0(dir_fun, "OneStep_H.R"))
	source(paste0(dir_fun, "weight_optimize.R"))
	source(paste0(dir_fun, "FlatTest_Hs.R"))
	```
 -->


## License

This project is licensed under the [MIT License](LICENSE). You are free to use, modify, and distribute the code with proper attribution.

<!-- ## Citation

If you find this work helpful, please cite the paper:


 -->

## Contact

For questions or feedback, please contact: \
Yutong Jin \
Fred Hutchinson Cancer Center\
yjin2@fredhutch.org\

