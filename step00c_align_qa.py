"""
step00c_align_qa.py

QA/QC for the Align CTD module, following Leah McRaven's slides, adapted to our instrument
(SBE 19plus V2). Align CTD shifts a channel in time relative to pressure so that temperature,
conductivity and oxygen describe the same parcel of water; getting it right minimises salinity
spiking and makes the down/up oxygen profiles agree.

Instrument-correct advances (from the training table, 19/19plus/19plus V2 row, NOT the 9plus row):
  - Temperature: +0.5 s relative to pressure  (applied as the fixed baseline here)
  - Oxygen (SBE 43 raw voltage): +3 to +7 s    (swept here to pick the best)
  - Conductivity: 0 on the 19plus (alignment is handled by the +0.5 s temperature advance);
    we sweep a small extra conductivity advance to confirm/tune it.

Metrics are measured in a WINDOW around the sharpest gradient (below the surface mixed layer),
because alignment spiking only shows at sharp interfaces; measuring over the whole cast lets real
near-surface down/up differences swamp the signal. Any "best" advance that lands on the edge of the
swept range is FLAGGED as boundary (not a real interior optimum), so a scattered or edge-pegged
result is not mistaken for a recommendation.

Per cast:
  1. Salinity assessment (conductivity timing): for each test conductivity advance, time-shift
     conductivity, recompute practical salinity (PSS-78), and measure salinity spiking as the std of
     the 24 Hz dS/dt inside the gradient window, plus the 2 db down/up salinity difference there.
  2. Oxygen assessment: for each test oxygen advance (3 to 7 s), time-shift SBE 43 voltage and
     measure down/up agreement vs potential temperature inside the window.
  3. Gradient strength: rank casts so tuning is done on the strongly-stratified ones.

The advance is emulated as a time-shift, so this needs no SBE binary and scans a whole folder. It
does not change any product: it tells you which advance to set in 02_alignctd.psa; the real Align
then runs in SBE.

USAGE: python step00c_align_qa.py
"""
from pathlib import Path
import re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---------------------------------------------------------------- paths / config
CTD_ROOT  = Path(r"C:\Projects\ctd_pipeline")
CRUISE    = "P45_06"
# 24 Hz cnvs with conductivity, temperature, pressure (and SBE43 voltage if present). Point this at
# your DatCnv (or Filter) output folder; the advance QC is done on 24 Hz data, before Align.
INPUT_DIR  = CTD_ROOT / "cruises" / CRUISE / "L1" / "CTD"
INPUT_GLOB = "*.cnv"
QC_PLOT_DIR = CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "qc_plots"; QC_PLOT_DIR.mkdir(parents=True, exist_ok=True)

FS = 24.0                       # scans per second
T_ADVANCE = 0.5                 # temperature advance (s), fixed, 19plus V2 instrument value
COND_ADVANCES = [-0.25, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.25]   # conductivity advances (s)
O2_ADVANCES   = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]               # oxygen advances (s)
PLOT_TOP_N = 3                  # detailed plots for the N strongest-gradient casts
BIN_DB = 2.0
SURF_SKIP_DB = 8.0              # ignore the top of the cast (mixed layer) when locating the gradient
GRAD_HALF_DB = 12.0             # half-width of the gradient window (centre +/- this), in db
MIN_GRAD_DSDP = 0.005           # below this 95th-pct |dS/dp| a cast is "weak gradient" (do not trust)

rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 9, "axes.titlesize": 10,
                 "axes.labelsize": 9, "axes.linewidth": 0.8, "xtick.direction": "in",
                 "ytick.direction": "in", "axes.grid": False, "legend.fontsize": 7.5,
                 "font.family": "DejaVu Sans"})

def _despine(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

# ---------------------------------------------------------------- seawater (PSS-78 + potential temp)
def sal78(R, T68, P):
    """PSS-78 practical salinity. R = conductivity ratio C/C(35,15,0), T68 in deg C (IPTS-68), P dbar."""
    a=[0.0080,-0.1692,25.3851,14.0941,-7.0261,2.7081]
    b=[0.0005,-0.0056,-0.0066,-0.0375,0.0636,-0.0144]
    c=[0.6766097,2.00564e-2,1.104259e-4,-6.9698e-7,1.0031e-9]
    d=[3.426e-2,4.464e-4,4.215e-1,-3.107e-3]
    e=[2.070e-5,-6.370e-10,3.989e-15]
    R=np.asarray(R,float)
    Rp=1+(P*(e[0]+e[1]*P+e[2]*P*P))/(1+d[0]*T68+d[1]*T68*T68+(d[2]+d[3]*T68)*R)
    rt=c[0]+c[1]*T68+c[2]*T68**2+c[3]*T68**3+c[4]*T68**4
    Rt=R/(Rp*rt); Rt=np.maximum(Rt,1e-12); s=np.sqrt(Rt)
    dS=(T68-15)/(1+0.0162*(T68-15))*(b[0]+b[1]*s+b[2]*Rt+b[3]*Rt*s+b[4]*Rt**2+b[5]*Rt**2*s)
    return a[0]+a[1]*s+a[2]*Rt+a[3]*Rt*s+a[4]*Rt**2+a[5]*Rt**2*s+dS

def practical_salinity(C_Sm, T90, P):
    """Practical salinity from conductivity [S/m], temperature [ITS-90 deg C], pressure [dbar]."""
    T68 = np.asarray(T90,float)*1.00024
    R = np.asarray(C_Sm,float) / 4.2914        # C(35,15,0) = 4.2914 S/m
    return sal78(R, T68, P)

def _adtg(S,T,P):
    a0,a1,a2,a3=3.5803e-5,8.5258e-6,-6.836e-8,6.6228e-10
    b0,b1=1.8932e-6,-4.2393e-8
    c0,c1,c2,c3=1.8741e-8,-6.7795e-10,8.733e-12,-5.4481e-14
    d0,d1=-1.1351e-10,2.7759e-12
    e0,e1,e2=-4.6206e-13,1.8676e-14,-2.1687e-16
    return (a0+(a1+(a2+a3*T)*T)*T + (b0+b1*T)*(S-35)
            + ((c0+(c1+(c2+c3*T)*T)*T)+(d0+d1*T)*(S-35))*P + (e0+(e1+e2*T)*T)*P*P)

def potential_temperature(S, T90, P, PR=0.0):
    """Potential temperature [deg C] referenced to PR dbar (UNESCO/Fofonoff RK4)."""
    S=np.asarray(S,float); T=np.asarray(T90,float)*1.00024; P=np.asarray(P,float)
    H=PR-P; XK=H*_adtg(S,T,P); T=T+0.5*XK; Q=XK; P2=P+0.5*H
    XK=H*_adtg(S,T,P2); T=T+0.29289322*(XK-Q); Q=0.58578644*XK+0.121320344*Q
    XK=H*_adtg(S,T,P2); T=T+1.707106781*(XK-Q); Q=3.414213562*XK-4.121320344*Q
    P3=P2+0.5*H; XK=H*_adtg(S,T,P3); T=T+(XK-2.0*Q)/6.0
    return T/1.00024                             # back to ITS-90

# ---------------------------------------------------------------- io / helpers
def read_cnv(path):
    txt=Path(path).read_text(encoding="latin-1",errors="replace").splitlines()
    names,ds={},None
    nr=re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]",re.I)
    for i,l in enumerate(txt):
        s=l.strip()
        if s.startswith("#") or s.startswith("*"):
            m=nr.search(l)
            if m: names[int(m.group(1))]=m.group(2).strip()
            if s.upper().startswith("*END*"): ds=i+1; break
    if ds is None: return None
    cols=[names[k] for k in sorted(names)]
    rows=[l.split()[:len(cols)] for l in txt[ds:] if l.strip() and len(l.split())>=len(cols)]
    return pd.DataFrame(rows,columns=cols).apply(pd.to_numeric,errors="coerce").reset_index(drop=True)

def pick(df,*cands):
    for c in cands:
        if c in df.columns: return c
    for c in df.columns:
        if any(x.lower() in c.lower() for x in cands): return c
    return None

def shift_advance(x, adv_s, fs=FS):
    """Advance a channel by adv_s seconds relative to pressure: value at scan i comes from a later
    scan (i + adv_s*fs), linearly interpolated. Edges hold the end values."""
    x=np.asarray(x,float); n=len(x); idx=np.arange(n)
    return np.interp(idx + adv_s*fs, idx, x, left=x[0], right=x[-1])

def descent_start(pv, win=int(FS), min_drop=0.10):
    pv=np.asarray(pv,float)
    if pv.size==0 or np.isnan(pv).all(): return 0,0
    imax=int(np.nanargmax(pv))
    if imax<win+1: return 0,imax
    seg=pv[:imax+1]; s=pd.Series(seg).rolling(win,center=True,min_periods=1).median().to_numpy()
    fwd=np.full(len(s),np.nan); fwd[:len(s)-win]=s[win:]-s[:len(s)-win]
    desc=np.where(np.isnan(fwd),True,fwd>min_drop); i=imax
    while i>0 and desc[i-1]: i-=1
    return i,imax

def bin_profile(p, v, step=BIN_DB):
    d=pd.DataFrame({"p":p,"v":v}).dropna()
    if d.empty: return np.array([]),np.array([])
    b=(d["p"]/step).round().astype(int)*step
    g=d.groupby(b)["v"].mean()
    return g.index.to_numpy(float), g.to_numpy(float)

def gradient_window(P, S, start, imax):
    """Locate the sharpest salinity gradient on the downcast (below the mixed layer) and return
    (lo, hi, centre, strength). strength = 95th pct |dS/dp| over the whole downcast. If the cast is
    too weakly stratified, centre is None and the window falls back to 'everything below SURF_SKIP'."""
    pb,sb=bin_profile(P[start:imax+1],S[start:imax+1])
    if len(pb)<5:
        return SURF_SKIP_DB, float(np.nanmax(P)), None, 0.0
    order=np.argsort(pb); pb,sb=pb[order],sb[order]
    dsdp=np.abs(np.gradient(sb,pb))
    strength=float(np.nanpercentile(dsdp,95))
    deep=pb>SURF_SKIP_DB
    if deep.sum()>=3 and strength>=MIN_GRAD_DSDP:
        centre=float(pb[deep][int(np.nanargmax(dsdp[deep]))])
        return centre-GRAD_HALF_DB, centre+GRAD_HALF_DB, centre, strength
    return SURF_SKIP_DB, float(np.nanmax(P)), None, strength

def _in(P, lo, hi):
    return (P>=lo) & (P<=hi)

# ---------------------------------------------------------------- assessments
def salinity_assessment(C,T,P, start,imax, lo,hi):
    """Per conductivity advance: dS/dt std (downcast, inside the gradient window) and 2 db down/up
    salinity RMS difference (inside the window). Returns (out, best, best_at_boundary)."""
    Ta=shift_advance(T,T_ADVANCE)
    dn=np.zeros(len(P),bool); dn[start:imax+1]=True
    up=np.zeros(len(P),bool); up[imax:]=True
    out={}
    for adv in COND_ADVANCES:
        Ca=shift_advance(C,adv); S=practical_salinity(Ca,Ta,P)
        m=dn & _in(P,lo,hi) & np.isfinite(S)
        dSdt=np.gradient(S[m]) if m.sum()>3 else np.array([np.nan])
        std=float(np.nanstd(dSdt))
        pbd,sbd=bin_profile(P[dn],S[dn]); pbu,sbu=bin_profile(P[up],S[up])
        common=np.intersect1d(pbd,pbu); common=common[(common>=lo)&(common<=hi)]
        if len(common):
            md={p:s for p,s in zip(pbd,sbd)}; mu={p:s for p,s in zip(pbu,sbu)}
            updown=float(np.sqrt(np.nanmean([(md[p]-mu[p])**2 for p in common])))
        else:
            updown=np.nan
        out[adv]={"S":S,"dSdt_std":std,"updown_rms":updown}
    best=min(out,key=lambda a: out[a]["dSdt_std"])
    at_bnd = best==min(COND_ADVANCES) or best==max(COND_ADVANCES)
    return out,best,at_bnd

def oxygen_assessment(O2v,S,T,P, start,imax, lo,hi):
    """Per oxygen advance: down/up agreement vs potential temperature inside the window (RMS of
    down-minus-up in theta bins). Returns (out, theta, best, best_at_boundary)."""
    theta=potential_temperature(S,T,P)
    dn=np.zeros(len(P),bool); dn[start:imax+1]=True
    up=np.zeros(len(P),bool); up[imax:]=True
    win=_in(P,lo,hi)
    out={}
    for adv in O2_ADVANCES:
        Oa=shift_advance(O2v,adv)
        thd,od=bin_profile(theta[dn&win],Oa[dn&win]); thu,ou=bin_profile(theta[up&win],Oa[up&win])
        common=np.intersect1d(np.round(thd,1),np.round(thu,1))
        if len(common):
            md={round(t,1):o for t,o in zip(thd,od)}; mu={round(t,1):o for t,o in zip(thu,ou)}
            rms=float(np.sqrt(np.nanmean([(md[t]-mu[t])**2 for t in common])))
        else:
            rms=np.nan
        out[adv]={"O2":Oa,"updown_rms":rms}
    valid={a:o for a,o in out.items() if np.isfinite(o["updown_rms"])}
    best=min(valid,key=lambda a: valid[a]["updown_rms"]) if valid else None
    at_bnd = best is not None and (best==min(O2_ADVANCES) or best==max(O2_ADVANCES))
    return out,theta,best,at_bnd

# ---------------------------------------------------------------- plots
def plot_salinity(cast,C,T,P,start,imax,sal_out,best,lo,hi,centre):
    Ta=shift_advance(T,T_ADVANCE); dn=slice(start,imax+1)
    fig,ax=plt.subplots(1,4,figsize=(12.5,6.0),sharey=True)
    cmap=plt.cm.viridis(np.linspace(0,1,len(COND_ADVANCES)))
    for c,adv in zip(cmap,COND_ADVANCES):
        S=sal_out[adv]["S"]; Sd=S[dn]; Pd=P[dn]
        ax[0].plot(Sd,Pd,lw=0.4,color=c,alpha=0.85,label=f"{adv:+.2f}s")
        ax[1].plot(np.gradient(Sd),Pd,lw=0.4,color=c,alpha=0.85)
        pb,sb=bin_profile(Pd,Sd); ax[2].plot(sb,pb,lw=1.0,color=c)
    Sb=sal_out[best]["S"]; pbd,sbd=bin_profile(P[dn],Sb[dn]); pbu,sbu=bin_profile(P[imax:],Sb[imax:])
    common=np.intersect1d(pbd,pbu)
    if len(common):
        md={p:s for p,s in zip(pbd,sbd)}; mu={p:s for p,s in zip(pbu,sbu)}
        ax[3].plot([md[p]-mu[p] for p in common],common,lw=0.9,color="#333333"); ax[3].axvline(0,color="k",lw=0.5,alpha=0.4)
    for a in ax:
        if centre is not None: a.axhspan(lo,hi,color="#F4A24A",alpha=0.12)   # gradient window
        a.invert_yaxis(); _despine(a)
    ax[0].set_ylabel("pressure (db)")
    ax[0].set_title("24 Hz salinity"); ax[0].set_xlabel("salinity (PSU)")
    ax[1].set_title("24 Hz dS/dt"); ax[1].set_xlabel("dS/dt (PSU/scan)")
    ax[2].set_title(f"{BIN_DB:g} db salinity"); ax[2].set_xlabel("salinity (PSU)")
    ax[3].set_title(f"{BIN_DB:g} db down-up (best {best:+.2f}s)"); ax[3].set_xlabel("\u0394 salinity (PSU)")
    ax[0].legend(title="cond advance",loc="lower left",fontsize=6.5,title_fontsize=7)
    win = "gradient window shaded" if centre is not None else "weak gradient: whole cast used"
    fig.suptitle(f"{cast}: Align CTD conductivity/salinity (T advance +{T_ADVANCE:g} s; {win})",fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.98])
    out=QC_PLOT_DIR/f"{cast}__qc_align_salinity.png"; fig.savefig(out,bbox_inches="tight"); plt.close(fig)
    return out

def plot_oxygen(cast,theta,P,start,imax,o2_out,best,lo,hi,centre):
    dn=slice(start,imax+1); up=slice(imax,len(P))
    fig,ax=plt.subplots(1,2,figsize=(9.0,6.0))
    cmap=plt.cm.plasma(np.linspace(0,0.9,len(O2_ADVANCES)))
    for c,adv in zip(cmap,O2_ADVANCES):
        O=o2_out[adv]["O2"]
        lbl=f"{adv:g}s" + (" *best" if adv==best else "")
        ax[0].plot(O[dn],P[dn],lw=0.5,color=c,label=lbl); ax[0].plot(O[up],P[up],lw=0.5,color=c,ls="--",alpha=0.7)
        ax[1].plot(O[dn],theta[dn],lw=0.5,color=c); ax[1].plot(O[up],theta[up],lw=0.5,color=c,ls="--",alpha=0.7)
    if centre is not None: ax[0].axhspan(lo,hi,color="#F4A24A",alpha=0.12)
    ax[0].invert_yaxis(); _despine(ax[0]); _despine(ax[1])
    ax[0].set_ylabel("pressure (db)"); ax[0].set_xlabel("SBE 43 voltage"); ax[0].set_title("O2 vs pressure (down solid, up dashed)")
    ax[1].set_ylabel("potential temperature (\u00b0C)"); ax[1].set_xlabel("SBE 43 voltage"); ax[1].set_title("O2 vs potential temperature")
    ax[0].legend(loc="best",fontsize=6.0,ncol=2)
    fig.suptitle(f"{cast}: Align CTD oxygen (down/up should agree near sharp changes)",fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.98])
    out=QC_PLOT_DIR/f"{cast}__qc_align_oxygen.png"; fig.savefig(out,bbox_inches="tight"); plt.close(fig)
    return out

# ---------------------------------------------------------------- run
def process_all(input_dir=INPUT_DIR, glob=INPUT_GLOB, plot_top_n=PLOT_TOP_N):
    files=sorted(Path(input_dir).glob(glob))
    if not files:
        print(f"no cnvs found in {input_dir} matching {glob}"); return
    rows=[]; cache={}
    for f in files:
        df=read_cnv(f)
        if df is None: continue
        pcol=pick(df,"prdM","prDM","prM"); tcol=pick(df,"tv290C","t090C","tv290")
        ccol=pick(df,"c0S/m","cond0S/m","c0mS/cm","c_S/m","cond")
        if pcol is None or tcol is None or ccol is None:
            print(f"{f.name}: missing P/T/C column, skipping"); continue
        ocol=pick(df,"sbeox0V","sbox0","oxygen raw","sbeoxV","sbeox0v")
        P=pd.to_numeric(df[pcol],errors="coerce").to_numpy(float)
        T=pd.to_numeric(df[tcol],errors="coerce").to_numpy(float)
        C=pd.to_numeric(df[ccol],errors="coerce").to_numpy(float)
        start,imax=descent_start(P)
        # window from the baseline salinity (T+0.5, C 0)
        S0=practical_salinity(shift_advance(C,0.0),shift_advance(T,T_ADVANCE),P)
        lo,hi,centre,strength=gradient_window(P,S0,start,imax)
        sal_out,best_c,cb=salinity_assessment(C,T,P,start,imax,lo,hi)
        row={"cast":f.stem,"n_scans":len(df),"grad_strength_dSdp":round(strength,4),
             "grad_centre_db":round(centre,1) if centre is not None else "",
             "gradient":"weak" if centre is None else "ok",
             "best_cond_advance_s":best_c,"cond_flag":"BOUNDARY" if cb else "interior",
             "dSdt_std_best":round(sal_out[best_c]["dSdt_std"],5),
             "cond_updown_rms":round(sal_out[best_c]["updown_rms"],4) if np.isfinite(sal_out[best_c]["updown_rms"]) else ""}
        o2_out=theta=best_o2=None; ob=False
        if ocol is not None:
            O2v=pd.to_numeric(df[ocol],errors="coerce").to_numpy(float)
            o2_out,theta,best_o2,ob=oxygen_assessment(O2v,S0,T,P,start,imax,lo,hi)
            row["best_o2_advance_s"]=best_o2 if best_o2 is not None else ""
            row["o2_flag"]="BOUNDARY" if ob else "interior"
            row["o2_updown_rms_best"]=round(o2_out[best_o2]["updown_rms"],4) if best_o2 is not None else ""
        rows.append(row)
        cache[f.stem]=dict(C=C,T=T,P=P,start=start,imax=imax,sal_out=sal_out,best_c=best_c,
                           o2_out=o2_out,theta=theta,best_o2=best_o2,has_o2=ocol is not None,
                           lo=lo,hi=hi,centre=centre)
    if not rows:
        print("no usable casts"); return
    summ=pd.DataFrame(rows).sort_values("grad_strength_dSdp",ascending=False).reset_index(drop=True)
    summ.to_csv(QC_PLOT_DIR/"qc_align_summary.csv",index=False,encoding="utf-8-sig")

    print(f"\nAlign CTD QC over {len(summ)} casts. Temperature advance fixed at +{T_ADVANCE:g} s (19plus V2).")
    print("Metrics measured in the gradient window; 'BOUNDARY' = best pegged at the sweep edge (not a")
    print("real optimum); 'weak' = too little stratification to judge. Trust interior results on")
    print("strong-gradient casts.\n")
    print(summ.to_string(index=False))

    # consensus among trustworthy casts (strong gradient + interior best)
    good=summ[(summ["gradient"]=="ok") & (summ["cond_flag"]=="interior")]
    if len(good):
        print(f"\nConductivity: interior (trustworthy) best advances -> {sorted(good['best_cond_advance_s'].tolist())}")
        print(f"  median = {np.median(good['best_cond_advance_s']):+.2f} s across {len(good)} casts")
    else:
        print("\nConductivity: no cast gave an interior optimum -> no evidence for a nonzero advance; keep 0 s.")
    if "o2_flag" in summ.columns:
        goodo=summ[(summ["gradient"]=="ok") & (summ["o2_flag"]=="interior")]
        if len(goodo):
            print(f"Oxygen: interior best advances -> {sorted(goodo['best_o2_advance_s'].tolist())}; "
                  f"median = {np.median(goodo['best_o2_advance_s']):.1f} s")
        else:
            print("Oxygen: no interior optimum -> judge from the O2-vs-potential-temperature plots; keep +5 s meanwhile.")

    made=[]
    for cast in summ["cast"].head(plot_top_n):
        c=cache[cast]
        made.append(plot_salinity(cast,c["C"],c["T"],c["P"],c["start"],c["imax"],c["sal_out"],c["best_c"],c["lo"],c["hi"],c["centre"]))
        if c["has_o2"] and c["o2_out"] is not None and c["best_o2"] is not None:
            made.append(plot_oxygen(cast,c["theta"],c["P"],c["start"],c["imax"],c["o2_out"],c["best_o2"],c["lo"],c["hi"],c["centre"]))
    print(f"\nsummary -> {QC_PLOT_DIR/'qc_align_summary.csv'}")
    print(f"{len(made)} plots (top {plot_top_n} gradient casts) -> {QC_PLOT_DIR}")
    print("Set the chosen values in 02_alignctd.psa (Conductivity, Oxygen raw); Temperature stays +0.5 s.")

if __name__ == "__main__":
    process_all()
