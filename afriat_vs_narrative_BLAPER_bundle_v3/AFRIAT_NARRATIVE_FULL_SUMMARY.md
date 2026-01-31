# Afriat Sets vs Narrative Rules — Full Summary (Artifacts Bundle)

Generated: 2026-01-31 06:12

## TL;DR one-liner
**Afriat gives a disciplined *set* of counterfactual behaviors; narrative/procedural rules give the *selection principle* that picks the point people actually do.**

---

## What we did
### Data / task
- Each subject faces repeated portfolio problems: allocate a fresh **$100** across three stocks (PapayaTech=X, AxolotlWorks=Y, QuipuQuantum=Z) under varying price vectors.
- History includes trials **t=1..49** (prices + dollar spends + cash leftover). The goal is to predict behavior at **t=50**.

### Two evaluation tracks
1) **Point-prediction track (t=1..49 → predict t=50):**
   - Baseline: mean-share scaling
   - Narrative: cross-validated selection over a library of interpretable rules
   - Afriat: concave utility with minimal additive slack (global δ), plus variants (δ_t, discretization, tie-break)

2) **Set-based Afriat track (t=1..25 → evaluate at t=50 as a set):**
   - Fit Afriat with **observation-specific slack δ_t** on early history (t=1..25).
   - Compute the Afriat counterfactual correspondence at t=50.
   - For visualization and diagnosis, project the set onto each stock’s spend coordinate and compute:
     - **Band**: [min spend, max spend]
     - **Coverage**: truth inside band
     - **Distance-to-band**: 0 if inside; else distance to nearest endpoint

---

## Main results — Point prediction (t=1..49 → t=50)
### Baseline vs narrative rule selection
| model                                       |   MAE_L1 |   Median_L1 |   P90_L1 |   Top1_acc |   NearAllIn_true>=90 |   NearAllIn_pred>=90 |
|:--------------------------------------------|---------:|------------:|---------:|-----------:|---------------------:|---------------------:|
| Baseline (mean-share scaling)               |    52.51 |       33.82 |   134.34 |       0.44 |                 0.09 |                 0    |
| Rule-based selection (CV over rule library) |    31.34 |       18.5  |    83.22 |       0.55 |                 0.09 |                 0.05 |

### Afriat variants (what helped, what didn’t)
| model                                                             |   MAE_L1 |   Median_L1 |   P90_L1 |   Top1_acc |   NearAllIn_pred>=90 |   NearAllIn_true>=90 |
|:------------------------------------------------------------------|---------:|------------:|---------:|-----------:|---------------------:|---------------------:|
| Rule-based (CV narrative rules)                                   |    31.34 |       18.5  |    83.22 |       0.55 |                 0.05 |                 0.09 |
| Afriat per-t δ_t + tie-break (closest to anchor) + 0.1-share grid |    47.8  |       31.82 |   115.88 |       0.54 |                 0.1  |                 0.09 |
| Baseline (mean-share scaling)                                     |    52.51 |       33.82 |   134.34 |       0.44 |                 0    |                 0.09 |
| Afriat per-t δ_t + 0.1-share postprocess                          |    52.66 |       43.67 |   111.56 |       0.49 |                 0.13 |                 0.09 |
| Afriat global δ + 0.1-share postprocess                           |    58.7  |       53.43 |   121.65 |       0.45 |                 0.18 |                 0.09 |
| Afriat + per-t δ_t (additive)                                     |    58.76 |       50.73 |   127.07 |       0.53 |                 0.08 |                 0.09 |
| Afriat + global δ (additive)                                      |    64.15 |       58.72 |   134.62 |       0.44 |                 0.14 |                 0.09 |

**Takeaway:** vanilla Afriat (global slack) is weak as a point predictor here; adding δ_t + tie-break + discretization helps, but rule-based narratives remain best on mean L1.

---

## Main results — Set evaluation (Afriat as a set, t=1..25 → t=50)
### Stock-by-stock Afriat band diagnostics at t=50
(Coverage and “narrative better|outside” are percentages.)
| stock        |   coverage |   median_width_$ |   median_dist_$ |   mean_dist_|outside_$ |   narrative_median_AE_$ |   narrative_better_|outside |
|:-------------|-----------:|-----------------:|----------------:|-----------------------:|------------------------:|----------------------------:|
| PapayaTech   |       23.4 |            4.695 |          11.145 |                22.2711 |                  11.485 |                        61   |
| AxolotlWorks |       20.1 |            4.42  |           9.53  |                20.5967 |                  10.575 |                        52   |
| QuipuQuantum |       18.8 |            4.39  |           8.02  |                18.5696 |                   9.495 |                        61.6 |

### Key X-only set fact (PapayaTech)
- Coverage (truth in band): **23.4%**
- Median band width: **$4.70**
- Median distance-to-band: **$11.14**
- Mean distance conditional on outside-band: **$22.27**
- Among outside-band cases, narrative point is closer than any Afriat-feasible X **61.0%** of the time.

---

## Narrative “types” of customers
### Selected using full history (t=1..49)
| family                         |   count |   pct |
|:-------------------------------|--------:|------:|
| Power-share (smooth)           |      55 |  35.7 |
| Case-based (kNN memory)        |      48 |  31.2 |
| Equal-dollar                   |      15 |   9.7 |
| Switching gate (ratio→kNN/reg) |      15 |   9.7 |
| Equal-share                    |      13 |   8.4 |
| Dominance (all-in)             |       5 |   3.2 |
| Share-anchor                   |       3 |   1.9 |

### Selected using early history only (t=1..25)
| family                  |   count |   pct |
|:------------------------|--------:|------:|
| Power-share (smooth)    |      61 |  39.6 |
| Case-based (kNN memory) |      41 |  26.6 |
| Equal-dollar            |      34 |  22.1 |
| Equal-share             |      13 |   8.4 |
| Dominance (all-in)      |       5 |   3.2 |

**Interpretation:** there is real procedural heterogeneity: some subjects follow smooth power-share behavior, others act case-by-case (memory), and a large minority use equal-dollar heuristics.

---

## What’s inside this bundle
### Papers
- `afriat_vs_narrative_paper.pdf` / `.tex` — original draft (X-focused)
- `afriat_vs_narrative_paper_v2.pdf` / `.tex` — updated draft with X, Y, Z set diagnostics + all figures

### Core prediction files (t=1..49 → t=50)
- `unified_predictions_all_subjects.json` — baseline mean-share scaling
- `unified_predictions_rule_based_selected.json` — narrative rule-selection predictions
- `unified_predictions_afriat_min_slack.json` — Afriat global slack predictions
- `rule_selection_by_subject.csv` — which narrative rule was selected per subject (full-history selection)

### Afriat set evaluation (t=1..25 → t=50)
- `afriat_band_papayatech_all_customers_t50.csv`
- `afriat_band_axolotlworks_all_customers_t50.csv`
- `afriat_band_quipuquantum_all_customers_t50.csv`
- `afriat_band_summary_all_stocks_t50.json` / `.csv`

### Plots (per stock)
For each stock (papayatech/axolotlworks/quipuquantum):
- `*_dist_to_band_hist_t50.png`
- `*_band_width_hist_t50.png`
- `*_dist_vs_ratio_t50.png`
- `*_width_vs_ratio_t50.png`

### One-customer deep dive
- `subject_920_*` and `Subject_920__*` files: per-period errors, set bands, and time-series figures.

### Source bundle
- `narrative_prompts_experiment_all_subjects_bundle.zip` (raw prompts/answers bundle used to generate these artifacts)

---

## Notes on reproducibility
- Some uploaded source files can expire in the environment over time. This bundle includes the source zip so you can re-run everything from the same inputs.
- All computations were deterministic given the bundle.