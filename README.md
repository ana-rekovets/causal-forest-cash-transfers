# Long-Term Impacts of Early-Life Cash Transfers: A Causal Forest Approach 

**Author:** Anastasiia Rekovets  

Replication and extension of heterogeneous treatment effect analysis of the early 20th-century **Mothers' Pension (MP) program** on child longevity. The analysis uses **Generalized Random Forests** (GRF) to estimate Conditional Average Treatment Effects (CATEs) and test for treatment effect heterogeneity across individual, family, and county characteristics.

> **Context:** The Mothers' Pension program was one of the first large-scale US welfare programs, providing cash transfers to single mothers. This analysis asks: did program acceptance increase children's lifespan, and for whom was the effect largest?

---

## Repository Structure

```text
.
├── scripts/
│   └── analysis_longevity_causal_forest.R   # Main analysis script
├── paper/
│   └── Rekovets_Mothers_Pension_Paper.pdf
├── presentation/
│   └── Rekovets_Mothers_Pension_Slides.pdf
├── README.md
└── data/                                    # Not included — see Data section
    ├── MP_1940.dta
    ├── MP_controls.dta
    ├── MP_data.dta
    ├── MP_Ohio.dta
    └── MP_WII.dta
```

---

## Methods

The script implements a two-stage **Causal Forest** (Wager & Athey, 2018) following the double-robustness / R-learner approach (Nie & Wager, 2021):

1. **Nuisance estimation** — separate regression forests for `E[Y|X]` and `E[W|X]`
2. **Variable selection** — initial forest used to identify covariates with above-average importance
3. **Tuned causal forest** — refitted on selected variables with `tune.parameters = "all"`
4. **Inference** — ATE, CATE, calibration tests, and doubly-robust heterogeneity tests

### Covariate Specifications

| Spec | Controls | County FE |
|------|----------|-----------|
| 1 | State + cohort FEs only | ✗ |
| 2 | All controls (continuous county) | ✗ |
| 3 | All controls | ✓ (main) |
| 4 | All controls + SSA outcome | ✓ |

### Robustness Checks

- Propensity score trimming (common support: 0.65–0.97)
- Forest without cluster adjustment
- Non-orthogonalized forest (`W.hat = mean(W)`)
- R-loss comparison across specifications
- County-level ANOVA on doubly-robust scores

---

## Requirements

| Package | Version |
|---------|---------|
| R | ≥ 4.0 |
| grf | ≥ 0.10.2 |
| haven | any |
| dplyr | any |
| ggplot2 | any |
| Hmisc | any |
| labelled | any |
| xtable | any |

Install all dependencies at once:

```r
install.packages(c("haven", "dplyr", "grf", "Hmisc",
                   "ggplot2", "labelled", "xtable"))
```

---

## Data

The raw `.dta` files are **not included** in this repository (proprietary / restricted access). Set the `base_path` variable at the top of the script to the folder containing the five Stata files:

```r
base_path <- "path/to/your/data/"
```

The dataset corresponds to the replication archive for:

> **[Original paper citation — add here]**  
> ICPSR Study 112988-V1

---

## Usage

```r
# 1. Set your data path in the script
# 2. Source the script
source("analysis_longevity_causal_forest.R")
```

All output (ATE tables, plots, LaTeX tables) is printed to the console or rendered inline. To save plots, wrap `print()` calls with `ggsave()` or redirect output to a PDF device.

---

## Key Results (example)

| Specification | ATE (log) | 95% CI | Effect (years) |
|---------------|-----------|--------|----------------|
| Spec 1: State + Cohort FEs | — | — | — |
| Spec 2: All controls, no county FE | — | — | — |
| Spec 3: All controls + county FE | — | — | — |
| Spec 4: SSA outcome | — | — | — |
| Trimmed sample | — | — | — |

*Fill in with your estimates after running the script.*

---

## References

- Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, 113(523), 1228–1242.
- Nie, X., & Wager, S. (2021). Quasi-oracle estimation of heterogeneous treatment effects. *Biometrika*, 108(2), 299–319.
- Tibshirani, J., Athey, S., Friedberg, R., Hadad, V., Hirshberg, D., Mayer, I., Sverdrup, E., Wager, S., & Wright, M. (2024). *grf: Generalized Random Forests* (R package).

---

## License

MIT License. See `LICENSE` for details.

