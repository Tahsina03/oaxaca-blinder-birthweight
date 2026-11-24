# Blinder-Oaxaca Decomposition of the Birth Weight Gap by Maternal Smoking

Babies of mothers who smoked during pregnancy weigh less at birth. This project uses a Blinder-Oaxaca decomposition to check how much of that gap comes from differences in the mothers' other characteristics, and how much is left over.

## Question

How much of the birth weight gap between non-smoking and smoking mothers is **explained** by observed characteristics, and how much is **unexplained**?

## Data

`birthwt` from the R `MASS` package: 189 births (115 non-smoking and 74 smoking mothers) at Baystate Medical Center, 1986. Mean birth weight was 3,056 g for non-smokers and 2,772 g for smokers.

## Methods

- Separate linear regressions of birth weight for non-smokers and smokers
- Threefold decomposition (endowments, coefficients, interaction)
- Twofold decomposition with pooled coefficients (Jann, 2008)
- Detailed decomposition by covariate
- Bootstrap confidence intervals (2,000 resamples, stratified by group)

## Results

| Component | Grams (95% CI) | % of gap |
|---|---|---:|
| **Total gap** | **283.8 (78.4, 476.8)** | 100 |
| Explained (pooled) | -38.7 (-184.4, 101.6) | -14 |
| Unexplained (pooled) | 322.5 (116.9, 521.1) | 114 |

![Threefold decomposition](output/figures/threefold.png)

- The observed characteristics do not explain the gap. Almost all of the 284 g gap is unexplained.
- Some risk factors (previous premature labour, uterine irritability, fewer early doctor visits) are more common in smokers. But the "other" race group, which has lower birth weights in this sample, is much more common in non-smokers (48% vs 16%). These work in opposite directions, so the explained part is close to zero.
- A large unexplained part fits with a direct effect of smoking, but it is not proof of causation. It also includes anything that was not measured. The sample is small, so the intervals are wide.

## How to run

Needs R 4.1+ with `MASS` and `ggplot2`.

```bash
Rscript R/analysis.R
```

Results are saved in `output/`. The seed is fixed.

## References

- Jann, B. (2008). The Blinder-Oaxaca decomposition for linear regression models. *The Stata Journal*, 8, 453–479.
- O'Donnell, O., van Doorslaer, E., Wagstaff, A., & Lindelow, M. (2008). *Analyzing Health Equity Using Household Survey Data*. World Bank.

## Author

Tahsina Fariha Ohi, MS Statistics, Ball State University
