# Tawam Cox / XGBoost-Cox Journal Reanalysis Package v1

## Purpose

This package reconstructs the Tawam survival-model analysis transparently because the original XGBoost-Cox script/output could not be recovered. It **does not attempt to reverse-engineer or tune a model to reproduce the old 0.850 number**.

It runs three prespecified models on the same 485-subject analytic cohort and the same frozen five-fold split (`MethodDev_CV_Fold5`, seed 20260808):

1. **Cox-5** — the verified Cox replication using five fitted terms: age 50–64, age >=65, CHD, diabetes, smoking.
2. **XGB-Cox-5 (matched-predictor sensitivity)** — XGBoost-Cox using exactly the same five terms as Cox. This isolates algorithm/functional-form differences.
3. **XGB-Cox-11 (journal reconstruction)** — XGBoost-Cox using the 11 baseline predictors stated in the manuscript: gender, continuous age, diabetes, CHD, vascular disease, smoking, hypertension, dyslipidemia, source obesity classification, ACEI/ARB exposure, and baseline eGFR.

The XGBoost configuration is fixed in the script and intentionally conservative for a cohort with only 56 events. It is a **post-hoc transparent reconstruction for the journal revision**, not the missing original model configuration, and it is not tuned against the held-out folds.

## Recommended system

- Windows 10/11
- Python 3.11 or 3.12 (64-bit)
- Internet access is needed only the first time, to install Python packages.

## Easiest way to run

1. Extract the ZIP to a normal local folder, for example:
   `D:\Hem\Tawam_Reanalysis_v1`
2. Double-click `RUN_TAWAM_ANALYSIS.bat`.
3. The script creates its own `.venv`, installs pinned dependencies, and runs 10,000 patient-bootstrap resamples.
4. When it says **COMPLETE**, open the `results` folder.
5. Zip the entire `results` folder and upload it back to ChatGPT.

If the full run is slow, double-click `RUN_TAWAM_ANALYSIS_FAST.bat` first. It uses 2,000 bootstrap resamples as a smoke test. For the manuscript, please run the full 10,000-resample version afterward.

## Outputs

The `results` folder will contain:

- `fold_metrics.csv` — fold-specific C-indices, mean and SD inputs.
- `overall_cindex_summary.csv` — pooled out-of-fold Harrell C with patient-bootstrap 95% CIs.
- `paired_cindex_differences.csv` — paired patient-bootstrap C-index differences.
- `oof_predictions.csv` — patient-level out-of-fold risk scores and 60-month predicted risks.
- `calibration_60m_summary.csv` — overall predicted vs Kaplan–Meier observed 60-month risk.
- `calibration_60m_by_quartile.csv` — risk-stratified calibration.
- `calibration_60m.png` — calibration plot.
- `cox_schoenfeld_trend_check.csv` — exploratory Schoenfeld-residual time-trend diagnostics for the five-term Cox model.
- `cox_schoenfeld_trend_check.png` — visual residual trends.
- `model_specification.json` — exact predictor lists and XGBoost hyperparameters.
- `run_metadata.json` — Python/package versions and input SHA-256.
- `RESULTS_SUMMARY.txt` — concise summary to send back to ChatGPT.

## Important interpretation rules

- Do **not** call XGB-Cox-11 a direct like-for-like comparison with Cox-5; it uses more predictors.
- Do **not** claim that the new XGBoost configuration is the original missing configuration.
- Do **not** replace manuscript values until the run completes and the outputs are reviewed.
- Calibration is descriptive because only 56 CKD progression events are available.
- The Schoenfeld residual checks are exploratory per-covariate time-trend diagnostics, not a formal global `cox.zph` test.

## Expected sanity checks

The input cohort should have:

- 485 subjects
- 56 events
- 429 right-censored observations
- five folds of 97 subjects each
- event counts 11, 11, 11, 11, 12

The verified Cox fold C-indices should be approximately:

`0.750, 0.803, 0.801, 0.795, 0.832`

Small floating-point differences are acceptable. If the Cox values differ materially, stop and send the results/error back before using the XGBoost outputs.
