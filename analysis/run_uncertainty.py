#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math
from datetime import date
from pathlib import Path
from statistics import NormalDist
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / 'config' / 'analysis_config.json').read_text())
DATA = json.loads((ROOT / 'inputs' / 'reported_summary.json').read_text())
OUT = ROOT / 'outputs'; FIG = ROOT / 'figures'
OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
seed = int(CFG['random_seed'])
rng = np.random.default_rng(seed)
alpha = 1 - float(CFG['confidence_level'])
z = NormalDist().inv_cdf(1-alpha/2)

def wilson(x:int,n:int):
    p=x/n
    d=1+z*z/n
    center=(p+z*z/(2*n))/d
    half=z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return max(0,center-half), min(1,center+half)

def metrics(tp,fp,fn,tn):
    tp=np.asarray(tp,dtype=float); fp=np.asarray(fp,dtype=float); fn=np.asarray(fn,dtype=float); tn=np.asarray(tn,dtype=float)
    precision=np.divide(tp,tp+fp,out=np.full_like(tp,np.nan),where=(tp+fp)>0)
    recall=np.divide(tp,tp+fn,out=np.full_like(tp,np.nan),where=(tp+fn)>0)
    accuracy=(tp+tn)/(tp+fp+fn+tn)
    specificity=np.divide(tn,tn+fp,out=np.full_like(tp,np.nan),where=(tn+fp)>0)
    f1=np.divide(2*precision*recall,precision+recall,out=np.full_like(tp,np.nan),where=(precision+recall)>0)
    return {'precision':precision,'recall':recall,'f1':f1,'accuracy':accuracy,'specificity':specificity}

c=DATA['social_test_confusion_matrix']
TP,FP,FN,TN=[int(c[k]) for k in ('TP','FP','FN','TN')]
N=TP+FP+FN+TN
point=metrics(TP,FP,FN,TN)
point={k:float(v) for k,v in point.items()}
point['prevalence']=(TP+FN)/N

B=int(CFG['bootstrap_replicates'])
p=np.array([TP,FP,FN,TN],dtype=float)/N
draws=rng.multinomial(N,p,size=B)
boot=metrics(draws[:,0],draws[:,1],draws[:,2],draws[:,3])
boot['prevalence']=(draws[:,0]+draws[:,2])/N

wilson_map={
 'precision':wilson(TP,TP+FP),
 'recall':wilson(TP,TP+FN),
 'accuracy':wilson(TP+TN,N),
 'specificity':wilson(TN,TN+FP),
 'prevalence':wilson(TP+FN,N)
}
rows=[]
for m in ['precision','recall','f1','accuracy','specificity','prevalence']:
    arr=np.asarray(boot[m],float)
    blo,bhi=np.nanquantile(arr,[alpha/2,1-alpha/2])
    w=wilson_map.get(m,(None,None))
    rows.append({'metric':m,'point_estimate':point[m],'wilson_95_low':w[0],'wilson_95_high':w[1],
      'parametric_bootstrap_95_low':float(blo),'parametric_bootstrap_95_high':float(bhi),
      'evidence_type':'empirical count reanalysis + simulation-based parametric bootstrap'})
with (OUT/'social_uncertainty.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

def d(s): return date.fromisoformat(s)
warning=DATA['server_warning_dates']
refs=[('principal_26_Feb_2020',CFG['principal_surveillance_milestone']),('alternative_30_Jan_2020',CFG['alternative_surveillance_milestone'])]
brows=[]
for rname,rdate in refs:
  for stage in ('early_pressure','confirmation'):
    delta=(d(rdate)-d(warning[stage])).days
    brows.append({'reference':rname,'reference_date':rdate,'stage':stage,'warning_date':warning[stage],
                  'lead_days_positive_before_reference':delta,
                  'interpretation':f'{delta} days before reference' if delta>=0 else f'{abs(delta)} days after reference',
                  'evidence_type':'exact reanalysis of archived dates'})
with (OUT/'benchmark_sensitivity.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=brows[0].keys()); w.writeheader(); w.writerows(brows)

ab=DATA['archived_pre_target_persistent_alert_days']; e=int(ab['early_pressure']); conf=int(ab['confirmation'])
ared={'early_pressure_days':e,'confirmation_days':conf,'absolute_reduction_days':e-conf,
      'relative_reduction':(e-conf)/e,'confirmation_to_early_ratio':conf/e,
      'evidence_type':'archived empirical server summary; descriptive only because days are serially dependent'}
with (OUT/'alert_burden_summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=ared.keys()); w.writeheader(); w.writerow(ared)

b=float(CFG['laplace_scale']); M=int(CFG['laplace_mc_replicates_per_margin'])
margins=[0.0,0.001,0.002,0.005,0.010,0.020]
prows=[]
for margin in margins:
    analytic=0.5*math.exp(-margin/b) if margin>0 else 0.5
    noise=rng.laplace(0,b,size=M)
    mc=float(np.mean(noise <= -margin))
    if margin <= 0.002: scenario='unfavorable_near_threshold'
    elif margin < 0.010: scenario='intermediate'
    else: scenario='favorable_well_separated'
    prows.append({'absolute_score_to_threshold_margin':margin,'scenario':scenario,
                  'analytic_single_release_flip_probability':analytic,
                  'monte_carlo_flip_probability':mc,'mc_replicates':M,
                  'evidence_type':'simulation/analytical sensitivity; not observed daily margins'})
with (OUT/'privacy_margin_sensitivity.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=prows[0].keys()); w.writeheader(); w.writerows(prows)

labels=['Precision','Recall','F1','Accuracy','Specificity']
sel={r['metric']:r for r in rows}
y=[sel[x.lower()]['point_estimate'] for x in labels]
lo=[y[i]-sel[labels[i].lower()]['parametric_bootstrap_95_low'] for i in range(len(labels))]
hi=[sel[labels[i].lower()]['parametric_bootstrap_95_high']-y[i] for i in range(len(labels))]
fig,ax=plt.subplots(figsize=(8,4.8))
ax.errorbar(range(len(labels)),y,yerr=[lo,hi],fmt='o',capsize=4)
ax.set_xticks(range(len(labels)),labels); ax.set_ylim(0.75,1.01); ax.set_ylabel('Metric value')
ax.set_title('Figure S1. Social-classifier uncertainty from parametric multinomial bootstrap')
ax.grid(axis='y',alpha=.25); fig.tight_layout(); fig.savefig(FIG/'figure_S1_social_uncertainty.png',dpi=220); plt.close(fig)

fig,ax=plt.subplots(figsize=(7.5,4.8))
xs=np.linspace(0,0.02,250); ys=0.5*np.exp(-xs/b)
ax.plot(xs,ys); ax.scatter([r['absolute_score_to_threshold_margin'] for r in prows],[r['analytic_single_release_flip_probability'] for r in prows])
ax.set_xlabel('Absolute score-to-threshold margin'); ax.set_ylabel('Single-release threshold-crossing probability')
ax.set_title('Figure S2. Laplace perturbation sensitivity at scale b = 0.002'); ax.set_yscale('log'); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig(FIG/'figure_S2_privacy_margin_sensitivity.png',dpi=220); plt.close(fig)

summary={'random_seed':seed,'social':{'N':N,'counts':{'TP':TP,'FP':FP,'FN':FN,'TN':TN},'point_estimates':point,
'parametric_bootstrap_replicates':B,'f1_bootstrap_95_CI':[sel['f1']['parametric_bootstrap_95_low'],sel['f1']['parametric_bootstrap_95_high']],
'precision_wilson_95_CI':[sel['precision']['wilson_95_low'],sel['precision']['wilson_95_high']],
'recall_wilson_95_CI':[sel['recall']['wilson_95_low'],sel['recall']['wilson_95_high']],
'accuracy_wilson_95_CI':[sel['accuracy']['wilson_95_low'],sel['accuracy']['wilson_95_high']]},
'benchmark_sensitivity':brows,'alert_burden':ared,'privacy_margin_sensitivity':prows,
'interpretation_limits':['The multinomial bootstrap quantifies sampling uncertainty conditional on the reported weak-label confusion matrix; it does not address weak-label construct validity.',
'Benchmark sensitivity changes the operational reference date, not the underlying warning trajectory.',
'Archived pre-target alert days are serially dependent and are therefore reported descriptively without an independence-based p-value.',
'Laplace margin analysis is a scenario/simulation sensitivity analysis; actual day-specific margins were unavailable and no empirical daily flip rate is inferred.']}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
