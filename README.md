### Parameter Estimation and Dynamic Prediction in Joint Models of Multiple Longitudinal Measures and Time-to-Event Outcome
A two-stage strategy is employed to jointly model multiple longitudinal measures and time-to-event outcomes. The initial stage entails estimating K one-marker joint models for each marker. In the subsequent stage, time-varying covariates are integrated into a proportional hazard model for parameter estimation. Both phases adhere to the Bayesian paradigm. Furthermore, this approach enables the computation of dynamic predictions based on the estimated parameter values.

### Installation
To acquire the latest development version of TSJM, you may utilize the following code snippet to install it directly from GitHub:

```
  # install.packages("devtools")
  devtools::install_github("tbaghfalaki/TSJM")
```
This will seamlessly fetch and install the most up-to-date version of TSJM for your use.


Note: Before installing this package, please install the *parallelsugar* package using the following command:

```
  # install.packages("devtools")
  devtools::install_github('nathanvan/parallelsugar')
```

### Example Usage
We consider two methods here. The first method is the full Bayesian approach, where both stages of the analysis are performed using Bayesian techniques. As a result, the second stage requires more computational time. Both stages can be easily executed using the TSJM package. Examples from the full Bayesian approach in TSJM can be found in the following:
- Parameter estimation: This analysis is presented [here](/Exam1.md)
- Dynamic prediction: This analysis is presented [here](/Exam2.md)

The second method uses the Rubin formula to estimate the standard deviation. We apply this approach to solve two examples: 
- Example 1: Without considring explanatroty variables  [here](/Exam3.md)
- Example 2: By considring explanatroty variables [here](/Exam4.md)


### Reference 
Baghfalaki, T., Hashemi, R., Helmer, C., & Jacqmin-Gadda, H. (2024). A two-stage joint modeling approach for multiple longitudinal markers and time-to-event data. arXiv preprint arXiv:2412.05765.
