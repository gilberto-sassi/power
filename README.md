# power

This function computes the power and sample size for basic testing hypothesis listed below.

## One population and variable

* Z test for mean (bilateral and unilateral): `pwr_z_test_1pop`
* t-test for mean (bilateral and unilateral): `pwr_t_test_1pop`
* Chi-squared for variance (bilateral and unilateral): `pwr_sigma_1pop`
* Proportion for test: `pwr_prop_1pop`

## Two populations and variables

* Z test for difference of means (bilateral and unilateral) with known variance: `pwr_z_test_2pop`
* F test to compare variance of two normal population: `pwr_sigma_2pop`
* t test for difference of means (bilateral and unilateral) with unknown variance and differences: `pwr_t_test_2pop_hetero`
* t test for difference of means (bilateral and unilateral) with unknown and equals: `pwr_t_test_2pop_homo`
* Paired t test: `pwr_paired_t_test`
* Test for proportions (bilateral and unilateral) in samples with more than 40 observations: `pwr_prop_2pop`

## Checking association

* Chi-sqaured to check association between two qualitative variables: `pwr_chisq_test_association` (implementation of power of test and sample size)
* Test for Pearson's correlation using Fisher's Z tranformation: `z_fisher_test` (implementation of test) and `pwr_z_fisher_test` (implementation of power of test and sample size)

## Anova

* Unbalanced ANOVA: `pwr_anova_balanced` (power of test and sample size)
* Balanced ANOVA: `pwr_anova_unbalanced` (only power of test)

## References

* MONTGOMERY, Douglas C.; RUNGER, George C. **Applied statistics and probability for engineers**. John Wiley & Sons, 2010.

