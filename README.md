# SW-MELTA: Sample-Weighted Mixture of Experts Latent Trait Analyzers

**SW-MELTA** extends the **Mixture of Experts Latent Trait Analyzers (MELTA)** framework to handle **complex survey designs**. MELTA is a model-based clustering approach for multivariate binary data that integrates finite mixtures with multidimensional latent traits. It further allows modeling the impact of covariates on cluster allocation, on the conditional responses’ distribution, on both (as in standard mixtures of experts models), or on neither. To deal with survey data, MELTA is extended by formally incorporating sampling weights into both the estimation algorithm and the uncertainty quantification. The proposed framework enables the identification of population-level latent clusters, while also assessing the role of observed features in cluster formation and observed responses in the target population.

## Repository Structure

* **`CODE.R`**: contains the main functions.
* **`main_sim.R`**: script to replicate the primary simulation study reported in the paper.
* **`additional_sim.R`**: script to replicate the supplementary simulation study.
* **`code_application.R`**: contains the main functions for the application.
* **`application.R`**: script to run model selection and to generate tables, plots, and final results for the empirical application.

* **`data.RData`**: dataset containing the binary indicators, covariates, and sampling weights.
* **`sw.melta.RData`**: main functions for the proposed **SW-MELTA** model and for the **MELTA** model.
* **`mlta.RData`**: main functions for the standard **MLTA** model.
* **`mlta.conc.RData`**: main functions for the **MLTA with concomitant variables** model.
