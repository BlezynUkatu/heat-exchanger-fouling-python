"""Illustrative counterflow exchanger; SI units. Run: python model.py."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent

def exchanger(rf=0., cold_flow=1.5, area=10., u_clean=500.):
    """Rf: combined area-normalized fouling resistance, m² K/W.
    Same reference area for U, Rf and A. Constant illustrative properties.
    """
    if rf < 0 or cold_flow <= 0 or area < 0 or u_clean <= 0:
        raise ValueError('Invalid flow, resistance, area or coefficient')
    ch, cc = 2.0*2200., cold_flow*4180.
    cmin, cmax = min(ch, cc), max(ch, cc)
    cr = cmin/cmax
    u = 1/(1/u_clean + rf)
    ntu = u*area/cmin
    eps = ntu/(1+ntu) if abs(cr-1)<1e-10 else -np.expm1(-ntu*(1-cr))/(1-cr*np.exp(-ntu*(1-cr)))
    q = eps*cmin*(120.-25.)
    return dict(rf_m2K_W=float(rf), cold_flow_kg_s=float(cold_flow), U_W_m2K=u,
                effectiveness=eps, Q_kW=q/1000., hot_out_C=120.-q/ch,
                cold_out_C=25.+q/cc)

def verify():
    """Independent LMTD root solve plus physical bounds and limiting cases."""
    errors=[]
    for rf in (0., .001, .002):
        for mc in (.5, 1.5, 3.):
            r=exchanger(rf, mc)
            ch,cc=4400.,mc*4180.
            def residual(q):
                d1,d2=95.-q/cc,95.-q/ch
                lm=(d1+d2)/2 if abs(d1-d2)<1e-8 else (d1-d2)/np.log(d1/d2)
                return q-r['U_W_m2K']*10*lm
            q_ref=brentq(residual,0.,min(ch,cc)*95*(1-1e-10))
            errors.append(abs(q_ref/1000-r['Q_kW']))
            assert np.isclose(q_ref/1000,r['Q_kW'],rtol=1e-9)
            assert 0<=r['effectiveness']<=1
            assert 25<=r['hot_out_C']<=120 and 25<=r['cold_out_C']<=120
            assert np.isclose(ch*(120-r['hot_out_C']),cc*(r['cold_out_C']-25))
    assert exchanger(area=0)['Q_kW']==0
    balanced=exchanger(cold_flow=4400/4180)
    assert np.isclose(balanced['effectiveness'],(5000/4400)/(1+5000/4400))
    assert exchanger(.001)['Q_kW']<exchanger(0)['Q_kW']
    return dict(independent_LMTD_cases=9, max_duty_difference_kW=max(errors),
                checks='Energy balance, temperature bounds, zero area, equal capacity rates, fouling trend: passed')

def main():
    for name in ('data','figures'): (ROOT/name).mkdir(exist_ok=True)
    rows=[exchanger(rf,mc) for rf in np.linspace(0,.002,51) for mc in np.linspace(.5,3,51)]
    np.savetxt(ROOT/'data/sensitivity.csv', [[r[k] for k in rows[0]] for r in rows],delimiter=',',header=','.join(rows[0]),comments='')
    clean,fouled=exchanger(),exchanger(.001)
    summary=dict(clean=clean,fouled=fouled,duty_reduction_percent=100*(1-fouled['Q_kW']/clean['Q_kW']),validation=verify(),cases=len(rows))
    (ROOT/'results.json').write_text(json.dumps(summary,indent=2))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'#f5f8fc','axes.facecolor':'#f5f8fc'})
    fig,ax=plt.subplots(figsize=(12,8))
    flows=np.linspace(.5,3,101)
    for rf,label,col in [(0,'Clean','#137e8b'),(.001,'Fouled: Rf = 0.001 m² K/W','#c25b35')]:
        ax.plot(flows,[exchanger(rf,m)['Q_kW'] for m in flows],label=label,color=col,lw=3)
    ax.set(xlabel='Cold-side mass flow (kg/s)',ylabel='Heat duty (kW)',title='How much heat can we recover?')
    ax.grid(alpha=.2); ax.legend(loc='lower right')
    fig.suptitle('PYTHON FOR CHEMICAL ENGINEERING  |  01',x=.12,ha='left',fontsize=13,color='#137e8b')
    fig.text(.12,.035,'Illustrative counterflow model • Constant clean U • No pressure-drop model\nHot inlet 120°C | Cold inlet 25°C | Area 10 m² | Blessing Ukatu',fontsize=10,color='#425466')
    fig.subplots_adjust(bottom=.18,top=.86)
    fig.savefig(ROOT/'figures/01_heat_recovery.png',dpi=180); plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,8))
    rs=np.linspace(0,.002,101)
    axs[0].plot(rs*1000,[exchanger(r)['Q_kW'] for r in rs],color='#137e8b',lw=3)
    axs[0].scatter([0,1],[clean['Q_kW'],fouled['Q_kW']],color='#c25b35',s=65,zorder=4)
    axs[0].set(xlabel='Fouling resistance (10⁻³ m² K/W)',ylabel='Heat duty (kW)',title='Less heat transferred')
    axs[1].plot(rs*1000,[exchanger(r)['hot_out_C'] for r in rs],label='Hot outlet',color='#c25b35',lw=3)
    axs[1].plot(rs*1000,[exchanger(r)['cold_out_C'] for r in rs],label='Cold outlet',color='#137e8b',lw=3)
    axs[1].set(xlabel='Fouling resistance (10⁻³ m² K/W)',ylabel='Outlet temperature (°C)',title='Outlet temperatures shift')
    for ax in axs: ax.grid(alpha=.2)
    axs[1].legend()
    fig.suptitle('FOULING: THE THERMAL PENALTY',fontsize=21,fontweight='bold',y=.96)
    fig.text(.08,.865,f"At Rf = 0.001 m² K/W: duty falls {summary['duty_reduction_percent']:.1f}% in this model",fontsize=16,color='#137e8b')
    fig.text(.08,.035,'Illustrative simulation, not plant measurements • Cold flow 1.5 kg/s • Fixed inlet conditions\nCombined fouling resistance on the same area basis as U • Blessing Ukatu',fontsize=10,color='#425466')
    fig.subplots_adjust(bottom=.18,top=.76,wspace=.32)
    fig.savefig(ROOT/'figures/02_fouling.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
