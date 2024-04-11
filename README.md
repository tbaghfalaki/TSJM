### Parameter Estimation and Dynamic Prediction in Joint Models of Multiple Longitudinal Measures and Time-to-Event Outcome
A two-stage strategy is employed to jointly model multiple longitudinal measures and time-to-event outcomes. The initial stage entails estimating K one-marker joint models for each marker. In the subsequent stage, time-varying covariates are integrated into a proportional hazard model for parameter estimation. Both phases adhere to the Bayesian paradigm. Furthermore, this approach enables the computation of dynamic predictions based on the estimated parameter values.

### Installation
To acquire the latest development version of TSJM, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/TSJM")
```
Note: Before installing this package, please install the *parallelsugar* package using the following command:

```
  # install.packages("devtools")
  devtools::install_github('nathanvan/parallelsugar')
```
This will seamlessly fetch and install the most up-to-date version of TSJM for your use.

### Example Usage

 > Examples from the TSJM can be found in the following:

- Parameter estimation: This analysis is presented [here](/Exam1.md)
- Dynamic prediction: This analysis is presented [here](/Exam2.md)

### Reference 
Baghfalaki, T., Hashemi, R., Helmer, C. & Jacqmin-Gadda, H. (2024). A Two-stage Joint Modeling Approach for Multiple Longitudinal Markers and Time-to-event Data. Submitted.
