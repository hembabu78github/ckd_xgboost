"""
Reproduce the Tawam method-development Cox analysis from tawam_analytic_cohort_v1.csv.
This script is intentionally independent of Excel parsing.
Dependencies: numpy, statsmodels, scikit-learn.
"""
import csv, math
import numpy as np
from sklearn.model_selection import StratifiedKFold
from statsmodels.duration.hazard_regression import PHReg

INPUT = "tawam_analytic_cohort_v1.csv"
SEED = 20260808

rows = list(csv.DictReader(open(INPUT, encoding="utf-8")))
time = np.array([float(r["TimeToEventMonths"]) for r in rows])
event = np.array([int(r["EventCKD35"]) for r in rows])
gender = np.array([int(r["GenderCode"]) for r in rows])
agecat = np.array([int(r["AgeCategoryCode"]) for r in rows])

X = np.column_stack([
    (agecat == 1).astype(int),
    (agecat == 2).astype(int),
    np.array([int(r["HistoryCHD"]) for r in rows]),
    np.array([int(r["HistoryDiabetes"]) for r in rows]),
    np.array([int(r["HistorySmoking"]) for r in rows]),
])
names = ["Age50_64", "Age65plus", "CHD", "Diabetes", "Smoking"]

def harrell_c(t, e, risk):
    c = 0.0
    ncomp = 0
    for i in range(len(t)):
        for j in range(i+1, len(t)):
            if e[i] == 1 and t[i] < t[j]:
                ncomp += 1
                c += 1 if risk[i] > risk[j] else 0.5 if risk[i] == risk[j] else 0
            elif e[j] == 1 and t[j] < t[i]:
                ncomp += 1
                c += 1 if risk[j] > risk[i] else 0.5 if risk[j] == risk[i] else 0
    return c/ncomp

fit = PHReg(time, X, status=event, ties="efron").fit()
print(fit.summary())
print("Apparent Harrell C:", harrell_c(time, event, X @ fit.params))

strata = event*2 + gender
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
cs = []
for fold, (tr, te) in enumerate(skf.split(np.zeros(len(rows)), strata), 1):
    m = PHReg(time[tr], X[tr], status=event[tr], ties="efron").fit(disp=0)
    c = harrell_c(time[te], event[te], X[te] @ m.params)
    cs.append(c)
    print(f"fold={fold} n={len(te)} events={event[te].sum()} c={c:.6f}")
print("Mean 5-fold C:", float(np.mean(cs)))
print("SD 5-fold C:", float(np.std(cs, ddof=1)))
