import csv,json,math
from pathlib import Path
R=Path(__file__).resolve().parents[1]
s=json.loads((R/'outputs'/'summary.json').read_text())
assert s['social']['N']==20000
assert s['social']['counts']=={'TP':346,'FP':3,'FN':68,'TN':19583}
assert abs(s['social']['point_estimates']['precision']-0.99140401146)<1e-9
assert abs(s['social']['point_estimates']['recall']-0.83574879227)<1e-9
assert s['benchmark_sensitivity'][0]['lead_days_positive_before_reference']==79
assert s['benchmark_sensitivity'][1]['lead_days_positive_before_reference']==24
assert s['alert_burden']['absolute_reduction_days']==55
assert abs(s['alert_burden']['relative_reduction']-55/65)<1e-12
print('Reproducibility checks passed.')
