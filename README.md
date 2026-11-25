# Reproducibility package: stage-specific federated decision fusion for outbreak warning

Public repository: https://github.com/elbitar24/outbreak-fusion-reproducibility

This package reproduces the **supplementary uncertainty and robustness analysis** added during final manuscript revision. It does not claim to reconstruct unavailable raw observations or rerun the original source-model training.

## Why uncertainty analysis was selected

An independent negative-control outbreak cohort was not available in the archived experimental package. Treating the 178 pre-event days inside the same COVID-19 episode as an independent negative-control cohort would overstate evidentiary independence. The lower-cost and more defensible reviewer response is therefore stronger uncertainty analysis using reported empirical counts, exact archived dates, archived alert-burden summaries, and clearly labeled simulations.

## One-command reproduction

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
./run_analysis.sh
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
.\run_analysis.ps1
```

The run regenerates:
- `outputs/social_uncertainty.csv`
- `outputs/benchmark_sensitivity.csv`
- `outputs/alert_burden_summary.csv`
- `outputs/privacy_margin_sensitivity.csv`
- `outputs/summary.json`
- `figures/figure_S1_social_uncertainty.png`
- `figures/figure_S2_privacy_margin_sensitivity.png`

## Statistical design

### Social classifier uncertainty
The reported confusion matrix (TP=346, FP=3, FN=68, TN=19,583; N=20,000) is reused without inventing post-level scores. Wilson 95% intervals are computed for binomially identifiable metrics. A 100,000-replicate parametric multinomial bootstrap (seed 20261125) propagates sampling uncertainty through precision, recall, F1, accuracy, specificity, and prevalence.

### Operational-reference sensitivity
Lead time is recomputed exactly for the principal 26 February 2020 CDC surveillance milestone and the alternative 30 January 2020 CDC person-to-person-spread milestone.

### Archived alert burden
The master-thesis server audit reported 65 pre-target persistent-alert days for early pressure and 10 for confirmation. The package reports the 55-day absolute and 84.6% relative reduction descriptively. No p-value is assigned because adjacent persistent-alert days are temporally dependent.

### Privacy-margin sensitivity
For Laplace scale b=0.002, the probability that a single perturbed scalar release crosses a threshold from an untreated margin m above the threshold is `0.5*exp(-m/b)`. Monte Carlo estimates validate the analytical calculation for unfavorable, intermediate, and favorable margins. Actual day-specific margins were unavailable, so this is a scenario analysis, not an empirical daily flip-rate estimate.

## Reproducibility boundaries
This package distinguishes original empirical results, deterministic reanalysis, and simulation-generated evidence. See `docs/PROVENANCE.md` and `docs/GAP_CLOSURE.md`.
