from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr
import sklearn
import statsmodels
from statsmodels.duration.hazard_regression import PHReg
import xgboost
from xgboost import XGBRegressor

SEED = 20260808
HORIZON_MONTHS = 60.0

BASE_DIR = Path(__file__).resolve().parent
DATA_PATH = BASE_DIR / "data" / "tawam_analytic_cohort_v1.csv"
RESULTS_DIR = BASE_DIR / "results"

COX5_NAMES = ["Age50_64", "Age65plus", "CHD", "Diabetes", "Smoking"]
XGB11_NAMES = [
    "GenderCode",
    "AgeBaselineYears",
    "HistoryDiabetes",
    "HistoryCHD",
    "HistoryVascular",
    "HistorySmoking",
    "HistoryHTN",
    "HistoryDLD",
    "HistoryObesity_Source",
    "ACEI_ARB",
    "eGFR_mL_min_1_73m2",
]

XGB_PARAMS = {
    "objective": "survival:cox",
    "eval_metric": "cox-nloglik",
    "n_estimators": 150,
    "max_depth": 2,
    "learning_rate": 0.03,
    "subsample": 0.80,
    "colsample_bytree": 0.80,
    "min_child_weight": 10,
    "reg_lambda": 5.0,
    "reg_alpha": 0.10,
    "random_state": SEED,
    "n_jobs": 1,
    "verbosity": 0,
    "tree_method": "hist",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def harrell_c(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    c = 0.0
    ncomp = 0
    n = len(time)
    for i in range(n):
        for j in range(i + 1, n):
            if event[i] == 1 and time[i] < time[j]:
                ncomp += 1
                c += 1.0 if risk[i] > risk[j] else 0.5 if risk[i] == risk[j] else 0.0
            elif event[j] == 1 and time[j] < time[i]:
                ncomp += 1
                c += 1.0 if risk[j] > risk[i] else 0.5 if risk[j] == risk[i] else 0.0
    return float(c / ncomp)


def comparable_pairs(time: np.ndarray, event: np.ndarray):
    early = []
    late = []
    for i in range(len(time)):
        for j in range(i + 1, len(time)):
            if event[i] == 1 and time[i] < time[j]:
                early.append(i); late.append(j)
            elif event[j] == 1 and time[j] < time[i]:
                early.append(j); late.append(i)
    return np.asarray(early, dtype=np.int32), np.asarray(late, dtype=np.int32)


def pair_scores(risk: np.ndarray, early: np.ndarray, late: np.ndarray) -> np.ndarray:
    a = risk[early]
    b = risk[late]
    return np.where(a > b, 1.0, np.where(a == b, 0.5, 0.0))


def bootstrap_cindices(
    time: np.ndarray,
    event: np.ndarray,
    risks: dict[str, np.ndarray],
    n_boot: int,
    seed: int = SEED,
):
    rng = np.random.default_rng(seed)
    n = len(time)
    early, late = comparable_pairs(time, event)
    scores = {k: pair_scores(v, early, late) for k, v in risks.items()}
    boots = {k: np.empty(n_boot, dtype=float) for k in risks}
    valid = 0
    while valid < n_boot:
        draw = rng.integers(0, n, size=n)
        counts = np.bincount(draw, minlength=n).astype(float)
        w = counts[early] * counts[late]
        den = w.sum()
        if den <= 0:
            continue
        for k in risks:
            boots[k][valid] = float(np.sum(w * scores[k]) / den)
        valid += 1
    return boots


def percentile_ci(x: np.ndarray, alpha=0.05):
    lo, hi = np.quantile(x, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)


def breslow_h0_at_horizon(time, event, hazard_ratio, horizon):
    h0 = 0.0
    event_times = np.sort(np.unique(time[(event == 1) & (time <= horizon)]))
    for tt in event_times:
        d = np.sum((time == tt) & (event == 1))
        denom = np.sum(hazard_ratio[time >= tt])
        if denom > 0:
            h0 += d / denom
    return float(h0)


def km_event_risk_at_horizon(time, event, horizon):
    surv = 1.0
    event_times = np.sort(np.unique(time[(event == 1) & (time <= horizon)]))
    for tt in event_times:
        at_risk = np.sum(time >= tt)
        d = np.sum((time == tt) & (event == 1))
        if at_risk > 0:
            surv *= (1.0 - d / at_risk)
    return float(1.0 - surv)


def fit_cox_fold(X, time, event, tr, te):
    fit = PHReg(time[tr], X[tr], status=event[tr], ties="efron").fit(disp=0)
    lp_tr = X[tr] @ fit.params
    lp_te = X[te] @ fit.params
    hr_te = np.exp(lp_te)
    h0_60 = float(np.asarray(fit.baseline_cumulative_hazard_function[0](HORIZON_MONTHS)))
    risk60 = 1.0 - np.exp(-h0_60 * hr_te)
    return fit, lp_te, risk60


def fit_xgb_fold(X, time, event, y_signed, tr, te):
    model = XGBRegressor(**XGB_PARAMS)
    model.fit(X[tr], y_signed[tr])
    hr_tr = model.predict(X[tr]).astype(float)
    hr_te = model.predict(X[te]).astype(float)
    h0_60 = breslow_h0_at_horizon(time[tr], event[tr], hr_tr, HORIZON_MONTHS)
    risk60 = 1.0 - np.exp(-h0_60 * hr_te)
    log_risk = np.log(np.maximum(hr_te, 1e-300))
    return model, log_risk, risk60


def make_calibration_table(df, risk60_dict, time, event):
    rows = []
    overall_obs = km_event_risk_at_horizon(time, event, HORIZON_MONTHS)
    summary = []
    for model, pred in risk60_dict.items():
        summary.append({
            "model": model,
            "n": len(pred),
            "mean_oof_predicted_60m_risk": float(np.mean(pred)),
            "km_observed_60m_risk": overall_obs,
            "absolute_gap": float(abs(np.mean(pred) - overall_obs)),
        })
        ranks = pd.Series(pred).rank(method="first")
        q = pd.qcut(ranks, 4, labels=False) + 1
        for quartile in range(1, 5):
            idx = np.where(q.to_numpy() == quartile)[0]
            rows.append({
                "model": model,
                "risk_quartile": quartile,
                "n": int(len(idx)),
                "events": int(event[idx].sum()),
                "mean_oof_predicted_60m_risk": float(np.mean(pred[idx])),
                "km_observed_60m_risk": km_event_risk_at_horizon(time[idx], event[idx], HORIZON_MONTHS),
            })
    return pd.DataFrame(summary), pd.DataFrame(rows)


def plot_calibration(cal_df: pd.DataFrame, path: Path):
    plt.figure(figsize=(6.8, 6.2))
    for model, g in cal_df.groupby("model", sort=False):
        plt.plot(g["mean_oof_predicted_60m_risk"], g["km_observed_60m_risk"], marker="o", label=model)
    lim = max(0.30, float(cal_df[["mean_oof_predicted_60m_risk", "km_observed_60m_risk"]].to_numpy().max()) * 1.08)
    plt.plot([0, lim], [0, lim], linestyle="--", linewidth=1, label="Ideal")
    plt.xlim(0, lim); plt.ylim(0, lim)
    plt.xlabel("Mean predicted 60-month CKD3–5 risk")
    plt.ylabel("Kaplan–Meier observed 60-month risk")
    plt.title("Out-of-fold 60-month calibration by risk quartile")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def schoenfeld_diagnostics(X5, time, event):
    fit = PHReg(time, X5, status=event, ties="efron").fit(disp=0)
    resid = np.asarray(fit.schoenfeld_residuals, dtype=float)
    event_idx = np.where(event == 1)[0]
    log_t = np.log(time[event_idx].astype(float))
    rows = []
    for j, name in enumerate(COX5_NAMES):
        r = resid[event_idx, j]
        mask = np.isfinite(r) & np.isfinite(log_t)
        rho, p = spearmanr(log_t[mask], r[mask])
        rows.append({
            "term": name,
            "n_events": int(mask.sum()),
            "spearman_rho_vs_log_event_time": float(rho),
            "p_value_exploratory": float(p),
            "note": "Exploratory Schoenfeld-residual time-trend diagnostic; not a formal global cox.zph test.",
        })
    return fit, pd.DataFrame(rows), resid


def plot_schoenfeld(schoen_df, path: Path):
    plt.figure(figsize=(8.0, 5.0))
    x = np.arange(len(schoen_df))
    plt.axhline(0, linewidth=1)
    plt.scatter(x, schoen_df["spearman_rho_vs_log_event_time"].to_numpy())
    plt.xticks(x, schoen_df["term"].tolist(), rotation=30, ha="right")
    plt.ylabel("Spearman rho: Schoenfeld residual vs log(event time)")
    plt.title("Cox proportional-hazards exploratory time-trend diagnostic")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bootstrap", type=int, default=10000, help="Number of patient bootstrap resamples (default 10000)")
    args = ap.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA_PATH)

    required = {
        "StudyID", "GenderCode", "AgeBaselineYears", "AgeCategoryCode",
        "HistoryDiabetes", "HistoryCHD", "HistoryVascular", "HistorySmoking",
        "HistoryHTN", "HistoryDLD", "HistoryObesity_Source", "ACEI_ARB",
        "eGFR_mL_min_1_73m2", "TimeToEventMonths", "EventCKD35", "MethodDev_CV_Fold5"
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing required columns: {sorted(missing)}")

    if len(df) != 485 or int(df["EventCKD35"].sum()) != 56:
        raise RuntimeError(f"Unexpected analytic cohort: n={len(df)}, events={int(df['EventCKD35'].sum())}; expected n=485/events=56")

    time = df["TimeToEventMonths"].to_numpy(float)
    event = df["EventCKD35"].to_numpy(int)
    fold = df["MethodDev_CV_Fold5"].to_numpy(int)
    agecat = df["AgeCategoryCode"].to_numpy(int)

    X5 = np.column_stack([
        (agecat == 1).astype(float),
        (agecat == 2).astype(float),
        df["HistoryCHD"].to_numpy(float),
        df["HistoryDiabetes"].to_numpy(float),
        df["HistorySmoking"].to_numpy(float),
    ])
    X11 = df[XGB11_NAMES].to_numpy(float)
    y_signed = np.where(event == 1, time, -time)

    expected_fold_sizes = [97] * 5
    fold_sizes = [int(np.sum(fold == k)) for k in range(1, 6)]
    fold_events = [int(event[fold == k].sum()) for k in range(1, 6)]
    if fold_sizes != expected_fold_sizes or fold_events != [11, 11, 11, 11, 12]:
        raise RuntimeError(f"Stored fold assignments failed sanity check: sizes={fold_sizes}, events={fold_events}")

    models = {
        "Cox-5": {"X": X5, "kind": "cox"},
        "XGB-Cox-5": {"X": X5, "kind": "xgb"},
        "XGB-Cox-11": {"X": X11, "kind": "xgb"},
    }

    oof_risk = {m: np.full(len(df), np.nan) for m in models}
    oof_risk60 = {m: np.full(len(df), np.nan) for m in models}
    fold_rows = []

    for k in range(1, 6):
        te = np.where(fold == k)[0]
        tr = np.where(fold != k)[0]
        for model_name, spec in models.items():
            X = spec["X"]
            if spec["kind"] == "cox":
                _, risk_te, risk60_te = fit_cox_fold(X, time, event, tr, te)
            else:
                _, risk_te, risk60_te = fit_xgb_fold(X, time, event, y_signed, tr, te)
            oof_risk[model_name][te] = risk_te
            oof_risk60[model_name][te] = risk60_te
            c = harrell_c(time[te], event[te], risk_te)
            fold_rows.append({
                "model": model_name,
                "fold": k,
                "n_test": int(len(te)),
                "events_test": int(event[te].sum()),
                "harrell_c": c,
            })

    fold_df = pd.DataFrame(fold_rows)
    fold_df.to_csv(RESULTS_DIR / "fold_metrics.csv", index=False)

    # Bootstrap pooled out-of-fold C-indices and paired differences.
    point_c = {m: harrell_c(time, event, r) for m, r in oof_risk.items()}
    boots = bootstrap_cindices(time, event, oof_risk, n_boot=args.bootstrap, seed=SEED)

    c_rows = []
    for m in models:
        lo, hi = percentile_ci(boots[m])
        fg = fold_df[fold_df.model == m]
        c_rows.append({
            "model": m,
            "fold_mean_c": float(fg.harrell_c.mean()),
            "fold_sd_c": float(fg.harrell_c.std(ddof=1)),
            "pooled_oof_c": point_c[m],
            "bootstrap_95ci_low": lo,
            "bootstrap_95ci_high": hi,
            "bootstrap_resamples": args.bootstrap,
        })
    c_df = pd.DataFrame(c_rows)
    c_df.to_csv(RESULTS_DIR / "overall_cindex_summary.csv", index=False)

    diff_rows = []
    pairs = [("XGB-Cox-11", "Cox-5"), ("XGB-Cox-5", "Cox-5"), ("XGB-Cox-11", "XGB-Cox-5")]
    for a, b in pairs:
        d = boots[a] - boots[b]
        lo, hi = percentile_ci(d)
        diff_rows.append({
            "contrast": f"{a} - {b}",
            "point_difference_pooled_oof_c": point_c[a] - point_c[b],
            "bootstrap_95ci_low": lo,
            "bootstrap_95ci_high": hi,
            "bootstrap_resamples": args.bootstrap,
        })
    diff_df = pd.DataFrame(diff_rows)
    diff_df.to_csv(RESULTS_DIR / "paired_cindex_differences.csv", index=False)

    # OOF predictions.
    oof = pd.DataFrame({
        "StudyID": df["StudyID"],
        "fold": fold,
        "TimeToEventMonths": time,
        "EventCKD35": event,
    })
    for m in models:
        safe = m.lower().replace("-", "_")
        oof[f"{safe}_oof_log_risk"] = oof_risk[m]
        oof[f"{safe}_oof_60m_risk"] = oof_risk60[m]
    oof.to_csv(RESULTS_DIR / "oof_predictions.csv", index=False)

    # Calibration at 60 months.
    cal_sum, cal_q = make_calibration_table(df, oof_risk60, time, event)
    cal_sum.to_csv(RESULTS_DIR / "calibration_60m_summary.csv", index=False)
    cal_q.to_csv(RESULTS_DIR / "calibration_60m_by_quartile.csv", index=False)
    plot_calibration(cal_q, RESULTS_DIR / "calibration_60m.png")

    # Full-cohort Cox coefficient table and Schoenfeld diagnostic.
    full_cox, schoen_df, _ = schoenfeld_diagnostics(X5, time, event)
    schoen_df.to_csv(RESULTS_DIR / "cox_schoenfeld_trend_check.csv", index=False)
    plot_schoenfeld(schoen_df, RESULTS_DIR / "cox_schoenfeld_trend_check.png")

    coef_df = pd.DataFrame({
        "term": COX5_NAMES,
        "coef_log_hazard": full_cox.params,
        "hazard_ratio": np.exp(full_cox.params),
        "se": full_cox.bse,
        "ci95_hr_low": np.exp(full_cox.params - 1.96 * full_cox.bse),
        "ci95_hr_high": np.exp(full_cox.params + 1.96 * full_cox.bse),
        "p_value": full_cox.pvalues,
    })
    coef_df.to_csv(RESULTS_DIR / "cox_full_cohort_coefficients.csv", index=False)

    model_spec = {
        "analysis_status": "post-hoc transparent journal reconstruction because original XGBoost-Cox code/output could not be recovered",
        "cv_split": {"source": "MethodDev_CV_Fold5", "seed_original": SEED, "fold_sizes": fold_sizes, "fold_event_counts": fold_events},
        "cox_5_predictors": COX5_NAMES,
        "xgb_cox_5_predictors": COX5_NAMES,
        "xgb_cox_11_predictors": XGB11_NAMES,
        "xgboost_fixed_parameters": XGB_PARAMS,
        "xgboost_selection_rule": "Fixed transparent reconstruction; no grid/random/Bayesian search and no optimization against held-out fold results.",
        "calibration_horizon_months": HORIZON_MONTHS,
    }
    (RESULTS_DIR / "model_specification.json").write_text(json.dumps(model_spec, indent=2), encoding="utf-8")

    metadata = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "statsmodels": statsmodels.__version__,
        "xgboost": xgboost.__version__,
        "matplotlib": matplotlib.__version__,
        "input_file": str(DATA_PATH.name),
        "input_sha256": sha256_file(DATA_PATH),
        "n_subjects": int(len(df)),
        "n_events": int(event.sum()),
        "n_censored": int((1-event).sum()),
        "bootstrap_resamples": args.bootstrap,
    }
    (RESULTS_DIR / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    lines = []
    lines.append("TAWAM JOURNAL REANALYSIS - RESULTS SUMMARY")
    lines.append("=" * 60)
    lines.append(f"Cohort: n={len(df)}, events={int(event.sum())}, censored={int((1-event).sum())}")
    lines.append(f"Stored fold event counts: {fold_events}")
    lines.append("")
    lines.append("Five-fold Harrell C (mean +/- SD) and pooled OOF C:")
    for _, r in c_df.iterrows():
        lines.append(
            f"  {r['model']}: fold mean={r['fold_mean_c']:.6f}, SD={r['fold_sd_c']:.6f}; "
            f"pooled OOF={r['pooled_oof_c']:.6f} (bootstrap 95% CI {r['bootstrap_95ci_low']:.6f}-{r['bootstrap_95ci_high']:.6f})"
        )
    lines.append("")
    lines.append("Paired pooled-OOF C-index differences:")
    for _, r in diff_df.iterrows():
        lines.append(
            f"  {r['contrast']}: {r['point_difference_pooled_oof_c']:.6f} "
            f"(95% CI {r['bootstrap_95ci_low']:.6f}-{r['bootstrap_95ci_high']:.6f})"
        )
    lines.append("")
    lines.append("60-month descriptive calibration:")
    for _, r in cal_sum.iterrows():
        lines.append(
            f"  {r['model']}: mean predicted={r['mean_oof_predicted_60m_risk']:.6f}; "
            f"KM observed={r['km_observed_60m_risk']:.6f}; absolute gap={r['absolute_gap']:.6f}"
        )
    lines.append("")
    lines.append("Cox PH exploratory Schoenfeld-residual time-trend diagnostics:")
    for _, r in schoen_df.iterrows():
        flag = "POTENTIAL TIME TREND" if r["p_value_exploratory"] < 0.05 else "no strong trend signal"
        lines.append(
            f"  {r['term']}: rho={r['spearman_rho_vs_log_event_time']:.4f}, p={r['p_value_exploratory']:.6g} [{flag}]"
        )
    lines.append("")
    lines.append("INTERPRETATION RULES:")
    lines.append("- XGB-Cox-11 is a post-hoc transparent reconstruction, not the missing original configuration.")
    lines.append("- XGB-Cox-11 is not a like-for-like comparison with Cox-5 because the predictor sets differ.")
    lines.append("- XGB-Cox-5 is the matched-predictor sensitivity analysis.")
    lines.append("- Calibration is descriptive because there are only 56 events.")
    lines.append("- Schoenfeld checks are exploratory per-term time-trend diagnostics, not a formal global cox.zph test.")
    (RESULTS_DIR / "RESULTS_SUMMARY.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print(f"\nOutputs written to: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
