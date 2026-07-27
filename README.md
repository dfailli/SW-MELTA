# SW-MELTA: Sample-Weighted Mixture of Experts Latent Trait Analyzers

**SW-MELTA** extends the **Mixture of Experts Latent Trait Analyzers (MELTA)** framework to handle **complex survey designs**. MELTA is a model-based clustering approach for multivariate binary data that integrates finite mixtures with multidimensional latent traits. It further allows modeling the impact of covariates on cluster allocation, on the conditional responses’ distribution, on both (as in standard mixtures of experts models), or on neither. To deal with survey data, MELTA is extended by formally incorporating sampling weights into both the estimation algorithm and the uncertainty quantification. The proposed framework enables the identification of population-level latent clusters, while also assessing the role of observed features in cluster formation and observed responses in the target population.

## Application

The empirical application focuses on measuring the **Digital Divide** across Italy using microdata from the 2021 edition of the *"Aspects of Daily Life"* (*Aspetti della vita quotidiana*) survey administered by **ISTAT** (Italian National Statistical Office). The aim is to identify population-level latent digital profiles in Italy based on multivariate binary indicators of digital competencies, while accounting for individuals' socio-demographic characteristics and survey design. 

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

## References

* **Failli, D., Marino, M.F., Martella, F.** (2025). *Mixture of experts latent trait analyzers*. In: D’Ambrosio, A., De Rooij, M., De Roover, K., Iorio, C., La Rocca, M. (eds) Supervised and Unsupervised Statistical Data Analysis. CLADAG-VOC 2025. Studies in Classification, Data Analysis, and Knowledge Organization. Springer, Cham. [https://doi.org/10.1007/978-3-032-03042-9_24](https://doi.org/10.1007/978-3-032-03042-9_24)
* **Failli, D., Marino, M.F., Martella, F.** (2024). *Finite mixtures of latent trait analyzers with concomitant variables for bipartite networks: An analysis of COVID-19 data*. Multivariate Behavioral Research, 59(4), 801–817.
* **Gollini, I., Murphy, T.B.** (2014). *Mixture of latent trait analyzers for model-based clustering of categorical data*. Statistics and Computing, 24(4), 569–588.
* **ISTAT** (2022). *Aspetti della vita quotidiana – Anno 2021. Nota metodologica*. Italian National Statistical Office. [Methodological Note PDF](https://www.istat.it/wp-content/themes/EGPbs5-child/microdata/download.php?id=/60/2021/01/Nota.pdf)

