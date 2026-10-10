"""Sharper rigorous eleven-term baseline; no lower-start splice is included.

Reuses the audited fixed-point/log and ordinary DDE enclosure backends. New:
exact logarithmic masses, cancellation b(R/a)/R=fbar(R/a)/a, concave
chords on one legal affine residual, and convex Fbar chords in U3/U4.
The frozen 184 table is used directly, never added as a second increment.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, time
from certify_conservative_chen import I,Q,log,display
from certify_extended_prefix import SieveTable
from certify_combined_chen_functional import SavedDefects,a,h,d,mu,e,ws,nu,C,boxes1

BASE=Path(__file__).resolve().parent
K=F(67619,10000)
OLDK=F(8459,1250)
def sha(name): return hashlib.sha256((BASE/name).read_bytes()).hexdigest()
def raw(n): return I(n,n,raw=True)

class Defects(SavedDefects):
    def __init__(self,name):
        self.raw=json.loads((BASE/name).read_text(encoding='utf-8'))
        assert int(self.raw['fixed_point_denominator'])==Q
        self.g=int(self.raw['grid_denominator'])
        assert self.g==50
        self.knots=[F(1)+F(i,self.g) for i in range(4*self.g+1)]
        self.H=list(map(int,self.raw['H_cell_lower_numerators']))
        self.h=list(map(int,self.raw['h_cell_lower_numerators']))
        assert len(self.H)==len(self.h)==4*self.g

class Ordinary:
    def __init__(self,g): self.t=SieveTable(g,8)
    def f(self,x):
        l,u=self.t.point(x.numerator,x.denominator)
        return I(l,u,raw=True)
    def upper(self,x):
        assert x>=1
        if x<=3:return I(1)/I(x)
        l,u=self.t.numerator(x.numerator,x.denominator,self.t.hl,self.t.hh,True)
        return I(l,u,raw=True)/I(x)
    def A(self,x): return I(x)*self.upper(x)
    def b(self,x): return I(x)*self.f(x)

def affine_range(p,t0,t1,u0,u1):
    c,b,d0=p
    return (c+b*(t0 if b>=0 else t1)+d0*(u0 if d0>=0 else u1),
            c+b*(t1 if b>=0 else t0)+d0*(u1 if d0>=0 else u0))

def pair_cell(t0,t1,u0,u1,enhance,diagonal,sf,df,cache):
    def ell(x,y):
        key=x,y
        if key not in cache:cache[key]=log(I(y/x))
        return cache[key]
    lt,lu=ell(t0,t1),ell(u0,u1)
    if diagonal:
        assert t0==u0 and t1==u1
        mass=lt*lt/2
        mt=I(t1)*lt-I(t1-t0)
        mu1=I(t1-t0)-I(t0)*lt
    else:
        mass=lt*lu;mt=I(t1-t0)*lu;mu1=I(u1-u0)*lt
    # Each tent is a LEGAL cap, minus the same t+u. Its active affine
    # plane may be used only if it is the smaller plane on the whole box.
    tents=[[(F(1,2),F(-1,2),F(-1)),(F(57,64),F(-5,2),F(-1))]]
    if enhance:
        tents.append([(F(1,2),F(-1,2),F(-1,2)),(F(57,64),F(-5,2),F(-5,2))])
        tents.append([(F(153,256)-F(3,8)*a,F(-1),F(-1))])
    lows=[];highs=[];affines=[]
    for tent in tents:
        ranges=[affine_range(p,t0,t1,u0,u1) for p in tent]
        lows.append(min(x[0] for x in ranges))
        highs.append(min(x[1] for x in ranges))
        if len(tent)==1:affines.extend(tent)
        else:
            for i in range(2):
                diff=tuple(tent[i][j]-tent[1-i][j] for j in range(3))
                if affine_range(diff,t0,t1,u0,u1)[1]<=0:affines.append(tent[i])
    rl,rh=max(lows),max(highs)
    assert rl>0
    best=(sf.f(rl/a)*mass*4/I(a)).lo
    # fbar is concave on [2,infinity). Its chord is a lower bound.
    # Use only one affine branch at a time; taking max of certified
    # box integrals is a choice of proof, not a sum of overlapping gains.
    for p in affines:
        pl,ph=affine_range(p,t0,t1,u0,u1)
        if pl<2*a or ph==pl:continue
        fl,fh=sf.f(pl/a),sf.f(ph/a)
        mom=I(p[0])*mass+I(p[1])*mt+I(p[2])*mu1
        val=(raw(fl.lo)*(I(ph)*mass-mom)+raw(fh.lo)*(mom-I(pl)*mass))/I(ph-pl)
        best=max(best,(val*4/I(a)).lo)
    hlow=df.lookup(I(rl/a,rh/a),True)
    defect=(4*hlow*mass/I(rh)).lo
    return max(0,best)+max(0,defect)

def upper_piece(l,r,A0,B0,n,sf,df):
    if l>=r:return I(0)
    # For s=max(2,R/a)<=3, A(s)=1 exactly, including the clamped part.
    cut=max(l,min(r,(A0-3*a)/B0));ordinary=I(0)
    if cut<r:
        ordinary+=4/I(A0)*(log(I(r/(A0-B0*r)))-log(I(cut/(A0-B0*cut))))
    if l<cut:
        step=(cut-l)/n
        for i in range(n):
            x=l+i*step;y=x+step
            sl=(A0-B0*y)/a;sr=(A0-B0*x)/a
            fl,fr=sf.upper(sl),sf.upper(sr)
            mass=log(I(y/x));moment=I(A0/a)*mass-I(B0/a*(y-x))
            # Convex Fbar lies below its endpoint chord.
            z=(raw(fl.hi)*(I(sr)*mass-moment)+raw(fr.hi)*(moment-I(sl)*mass))/I(sr-sl)
            ordinary+=4*z/I(a)
    return ordinary-df.correction(l,r,A0,B0,True)

def uniform_four(K0):
    L=log(I(h/a));J=[]
    series=I(1);factorial=1
    for n in range(4):
        if n:
            factorial*=n;series+=L**n/I(factorial)
        J.append(I(factorial)*(1-I(a/h)*series))
    u10=I(K0*C/(2*a))*(L**2*J[1]-2*L*J[2]+J[3])
    u11=I(K0*C/a)*log(I(d/h))*(L*J[1]-J[2])
    return u10,u11

def run(n,lines,g,bank,outname):
    started=time.monotonic();sf=Ordinary(g);df=Defects(bank);terms={}
    assert a>=F(25,352) and F(1,12)<a<h<mu<d<e and h<2*a
    assert C*(1/h-4)>1 and C*((1-d)/h-3)>1
    terms['L1']=4*sf.b(nu/a)/I(nu)
    sstar=F(249,50);sh=(1+h)/(2*h)
    terms['L2']=8*sf.b(sh)-4*sf.b(sstar)/I(sstar)+df.l2_correction()
    first=upper_piece(a,mu,F(1,2),F(1,2),lines,sf,df)
    terms['U3']=first+upper_piece(mu,e,F(57,64),F(5,2),lines,sf,df)+upper_piece(e,F(1,3),F(1,2),F(1),lines,sf,df)
    terms['U4']=first+upper_piece(mu,d,F(57,64),F(5,2),lines,sf,df)
    print(json.dumps({'stage':'ordinary_and_one_prefix','seconds':time.monotonic()-started}),flush=True)
    step=(h-a)/n;cache={};v=[a+i*step for i in range(n+1)]
    total=0
    for i in range(n):
        for j in range(i+1):total+=pair_cell(v[i],v[i+1],v[j],v[j+1],True,i==j,sf,df,cache)
    terms['L5']=raw(total)
    print(json.dumps({'stage':'L5','lower':total/Q,'seconds':time.monotonic()-started}),flush=True)
    total=0
    for left,right,enhance in [(h,mu,True),(mu,d,False)]:
        dt=(right-left)/n
        for i in range(n):
            t0=left+i*dt
            for j in range(n):total+=pair_cell(t0,t0+dt,v[j],v[j+1],enhance,False,sf,df,cache)
        print(json.dumps({'stage':'L6_piece','enhance':enhance,'seconds':time.monotonic()-started}),flush=True)
    terms['L6']=raw(total)
    def ainterval(t):
        lo=F(Q,t.hi)-1;hi=F(Q,t.lo)-1
        return I(sf.A(lo).lo,sf.A(hi).hi,raw=True)
    terms['U7']=K*boxes1(lambda t:(log((1-2*t)/t)+ainterval(t)-1)/(t*(1-t)),d,F(1,3),lines*8)
    terms['U8']=K*boxes1(lambda t:log(2-3*t)/(t*(1-t)),a,F(1,3),lines*8)
    terms['U9']=K*boxes1(lambda t:log((1-t-d)/I(d))/(t*(1-t)),h,d,lines*8)
    # To retain the old patch without changing its reference measure, use
    # the old outer K only in these two uniform terms. A later new patch
    # can instead use uniform_four(K) and replace the old patch entirely.
    terms['U10'],terms['U11full']=uniform_four(OLDK)
    c=(3*terms['L1']+terms['L2']+terms['L5']+terms['L6']-terms['U3']-terms['U4']-2*terms['U7']-terms['U8']-terms['U9']-terms['U10']-terms['U11full'])/4
    patch=json.loads((BASE/'combined-four-prefix-min-n32.json').read_text(encoding='utf-8'))
    oldgain=F(patch['Chen_gain_lower']);final=F(c.lo,Q)+oldgain
    previous=F(json.loads((BASE/'certify_variable_kernel_global_increment.json').read_text(encoding='utf-8'))['raw_c_lower'])
    out={'status':'EXACT_SHARPER_FULL_ELEVEN_TERM_BASELINE_UNDER_EXISTING_INSERTION_INPUTS','parameters':{'a':str(a),'h':str(h),'d':str(d),'K_U7_U9':str(K),'K_uniform_U10_U11':str(OLDK)},'bank':bank,'meshes':{'pairs':n,'ordinary':g,'one_prefix':lines,'switch_line':lines*8},'Q':str(Q),'terms':{k:display(v) for k,v in terms.items()},'uniform_c_lower':str(F(c.lo,Q)),'preserved_old_patch_once':str(oldgain),'raw_c_lower':str(final),'raw_c_diagnostic':float(final),'old_184_c_lower':str(previous),'improvement_over_old_184_lower':str(final-previous),'improvement_diagnostic':float(final-previous),'scope':['Complete old eleven-term baseline only; no lower-start three/four/five-prefix increment.','Same frozen bank used directly; no duplicate table delta.','One affine branch per pair box, concave fbar chord, exact logarithmic rectangle or diagonal triangle mass.','Convex Fbar chords in upper one-prefix integrals; exact clamped A=1 region.','Old four-prefix patch kept once against its original outer K; all other K use 6.7619.','Shared rigorous I/log and SieveTable primitives; analytic insertion/distribution hypotheses remain inherited.','Positive-term fields report proved lower bounds, not upper enclosures of the actual table functional.'],'elapsed_seconds':time.monotonic()-started,'bindings_sha256':{x:sha(x) for x in [Path(__file__).name,'certify_conservative_chen.py','certify_extended_prefix.py','certify_combined_chen_functional.py',bank,'combined-four-prefix-min-n32.json','certify_variable_kernel_global_increment.json','finite-table-chen-insertion.tex']}}
    (BASE/outname).write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'result':outname,'c':float(final),'gain':float(final-previous),'seconds':out['elapsed_seconds']}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pairs',type=int,default=128);ap.add_argument('--lines',type=int,default=2048);ap.add_argument('--ordinary',type=int,default=4096);ap.add_argument('--bank',default='certify_variable_kernel_finite_table.json');ap.add_argument('--output',default='sharp-184-baseline-128.json');x=ap.parse_args();run(x.pairs,x.lines,x.ordinary,x.bank,x.output)
