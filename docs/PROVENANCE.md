# Provenance and evidence classes

The package deliberately distinguishes three evidence classes.

1. **Original empirical/archived results**: social test confusion counts (N=20,000; TP=346, FP=3, FN=68, TN=19,583), archived warning dates (9 Dec 2019 and 2 Feb 2020), privacy-utility summary (MAD 0.001754; maximum 0.009493; correlation 0.999958), and archived pre-target persistent-alert-day burden (65 early-pressure days, 10 confirmation days).
2. **New deterministic statistical reanalysis**: Wilson 95% confidence intervals, date arithmetic under the alternative 30 Jan 2020 operational milestone, and alert-burden reduction (55 days; 84.6%).
3. **New simulation-based sensitivity analysis**: a parametric multinomial bootstrap conditional on the reported social confusion matrix, and Laplace threshold-crossing probabilities conditional on hypothetical score-to-threshold margins. These simulations do not create or replace empirical observations.

Source documents used for provenance:
- `Paper_Revised_Final_High_Level_Redlined.docx` (latest revised manuscript available during analysis).
- `CENG695 - Master_Thesis Part B - Malek Srouji and Mohamad Hajj.pdf`, Chapter 7, especially Tables 7-9 to 7-11 and the critical evaluation (printed pp. 77-79; PDF pages 88-90).
