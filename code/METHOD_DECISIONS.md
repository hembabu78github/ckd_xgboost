# Method decisions for the journal reconstruction

## Evidence boundary

The original XGBoost-Cox code/output used to report a mean five-fold C-index around 0.850 could not be recovered. The reconstruction therefore starts from the frozen public-data analytic cohort and stored CV fold assignments. It does not search for hyperparameters that recreate the old number.

## Cox model

The Cox model reproduces the existing method-development analysis with five fitted terms:

- Age 50–64 indicator
- Age >=65 indicator
- History of coronary heart disease
- History of diabetes
- History of smoking

This differs from the 11-predictor XGBoost model described in the manuscript. The new package makes this difference explicit.

## XGBoost-Cox models

Two XGBoost-Cox analyses are run:

- **XGB-Cox-5:** same five predictors as Cox, for a matched-predictor sensitivity analysis.
- **XGB-Cox-11:** the 11 baseline predictors previously stated in the manuscript, as a transparent reconstructed exploratory model.

The fixed XGBoost configuration is deliberately shallow and regularized. No grid search, random search, Bayesian optimization, or held-out-fold selection is performed.

## C-index uncertainty

Out-of-fold patient predictions are pooled, and 95% CIs are estimated by patient bootstrap. Paired differences are calculated from the same bootstrap samples so within-patient pairing is preserved.

## Calibration

For each training fold, a baseline cumulative hazard is derived from the training data. A 60-month event probability is then calculated for the held-out fold. Calibration is summarized overall and by out-of-fold risk quartile using Kaplan–Meier observed 60-month risk. This is descriptive because the number of events is small.

## Proportional-hazards diagnostic

The full-cohort five-term Cox model produces Schoenfeld residuals. Spearman correlation between each event-time residual series and log(event time) is reported as an exploratory time-trend diagnostic. This is not represented as a formal global Grambsch–Therneau/cox.zph test.
