"""Combined 56-generation W1/W2 table in the eleven-term functional.

Adapted from the independently audited independent_finite_chen_functional.py. Reuse its ordinary DDE,
closed-cell lookup, physical pair domains, exact one-prefix corrections,
and all five switched kernels. This is an arithmetic certificate conditional
on the separately audited insertion/distribution/seed interfaces.
"""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left, bisect_right
import argparse, hashlib, json, time
from certify_conservative_chen import I,Q,log,display

BASE=Path(__file__).resolve().parent
a=F(417,5000); h=F(19,125); d=F(59,250); K=F(8459,1250)
mu=F(25,128); e=F(57,224); ws=F(25,224); nu=F(249,448)
C=F(70893,125000)

def sha(p):return hashlib.sha256((BASE/p).read_bytes()).hexdigest()
def low_min(x,y):return I(min(x.lo,y.lo),min(x.hi,y.hi),raw=True)
def cap(t):return low_min((1+t)/2,I(F(57,64))-F(3,2)*t)
def boxes1(f,l,r,n):
    step=(r-l)/n;lower=upper=0
    for i in range(n):
        z=f(I(l+i*step,l+(i+1)*step));lower+=z.lo;upper+=z.hi
    return I(lower,upper,raw=True)*I(step)

class Delay:
    def __init__(self,g):
        self.g=g;self.A=[I(1)];self.b=[I(0)]
        for k in range(1,6*g+1):
            right=F(g+k,g);left=right-F(1,g)
            if right<=3:A=I(1)
            else:
                lag=I(self.b[k-g-1].lo,self.b[k-g].hi,raw=True)
                A=self.A[-1]+lag*log(I((right-1)/(left-1)))
            if right<=2:b=I(0)
            elif right<=4:b=log(I(right-1))
            else:
                lag=I(self.A[k-g-1].lo,self.A[k-g].hi,raw=True)
                b=self.b[-1]+lag*log(I((right-1)/(left-1)))
            self.A.append(A);self.b.append(b)
        assert all(x.lo>=0 and x.lo<=x.hi for x in self.A+self.b)
    def at(self,x,which='A'):
        assert x.hi<=7*Q
        if which=='b' and x.hi<=2*Q:return I(0)
        if which=='A' and x.hi<=3*Q:return I(1)
        arr=self.A if which=='A' else self.b
        lower=max(0,(x.lo-Q)*self.g//Q)
        upper=-(-((x.hi-Q)*self.g)//Q)
        assert upper<len(arr)
        return I(arr[lower].lo,arr[upper].hi,raw=True)

class SavedDefects:
    def __init__(self):
        self.raw=json.loads((BASE/'certify_combined_finite_table.json').read_text())
        assert self.raw['status']=='EXACT_FINITE_COMBINED_W1_W2_BUCHSTAB_ARITHMETIC_UNDER_ACTION_AND_TYPED_INPUTS'
        assert int(self.raw['fixed_point_denominator'])==Q
        for n,digest in self.raw['dependencies_sha256'].items():assert sha(n)==digest
        assert sha('certify_combined_finite_table.py')==self.raw['source_sha256']
        assert self.raw['total_finite_generations']==56
        audit=json.loads((BASE/'combined-finite-table-independent-audit.json').read_text())
        assert audit['status']=='PASS_INDEPENDENT_ALL_GENERATION_FULL_ARRAY_COMBINED_RECURRENCE_REPLAY'
        for name,digest in audit['bindings_sha256'].items():assert sha(name)==digest,name
        self.g=int(self.raw['grid_denominator']);self.knots=[F(1)+F(i,self.g) for i in range(4*self.g+1)]
        self.H=list(map(int,self.raw['H_cell_lower_numerators']))
        self.h=list(map(int,self.raw['h_cell_lower_numerators']))
    def lookup(self,x,lower=False):
        if x.lo<Q or x.hi>5*Q or (lower and x.lo<2*Q):return I(0)
        # Closed cells meeting [x.lo/Q,x.hi/Q], including both sides at a knot.
        start=max(0,-(-((x.lo-Q)*self.g)//Q)-1)
        stop=min(4*self.g-1,(x.hi-Q)*self.g//Q)
        nums=self.h if lower else self.H
        return I(min(nums[start:stop+1]),min(nums[start:stop+1]),raw=True)
    def correction(self,l,r,A0,B0,clamped=False):
        """4 int H(max(2,(A0-B0*t)/a))/(t*(A0-B0*t)) dt."""
        total=I(0)
        def mass(x,y):
            return (log(I(y/(A0-B0*y)))-log(I(x/(A0-B0*x))))/I(A0)
        if clamped:
            cl=max(l,(A0-2*a)/B0)
            if cl<r:total+=4*self.lookup(I(2))*mass(cl,r)
            r=min(r,(A0-2*a)/B0)
        for i,num in enumerate(self.H):
            if num==0:continue
            left=max(l,(A0-a*self.knots[i+1])/B0)
            right=min(r,(A0-a*self.knots[i])/B0)
            if left<right:total+=I(num,num,raw=True)*4*mass(left,right)
        return total
    def l2_correction(self):
        total=I(0)
        for i,num in enumerate(self.H):
            left=max(ws,1/(2*self.knots[i+1]+1));right=min(h,1/(2*self.knots[i]+1))
            if left<right:
                weight=8*(log(I(right/(1-right)))-log(I(left/(1-left))))
                total+=I(num,num,raw=True)*weight
        return total

def run(pair_n=384,line_n=8192,g=8192):
    start=time.monotonic()
    assert F(1,12)<a<h<mu<d<e and ws<h<2*a
    assert F(1,3)<h+d<F(1,2) and 4*h+d<1 and 4*a+d<1
    assert a>=F(25,352) and d>=F(1,5)
    assert C*(1/h-4)>1 and C*((1-d)/h-3)>1
    saved=json.loads((BASE/'certify_taylor_fifth_seed_transfer.json').read_text())
    assert K==F(saved['safe_K'])
    for n,digest in saved['dependencies_sha256'].items():assert sha(n)==digest
    sf=Delay(g);df=SavedDefects();terms={}
    terms['L1']=4*sf.at(I(nu/a),'b')/I(nu)
    sstar=F(249,50);sh=(1+h)/(2*h)
    terms['L2']=8*sf.at(I(sh),'b')-4*sf.at(I(sstar),'b')/I(sstar)+df.l2_correction()
    print(json.dumps({'stage':'delay_and_exact_L2','seconds':time.monotonic()-start}),flush=True)
    def upper(l,r,A0,B0):
        if l>=r:return I(0)
        ordinary=boxes1(lambda t:4*sf.at(((A0-B0*t)/I(a)).maximum(2))/(t*(A0-B0*t)),l,r,line_n)
        return ordinary-df.correction(l,r,A0,B0,True)
    terms['U3']=upper(a,mu,F(1,2),F(1,2))+upper(mu,e,F(57,64),F(5,2))+upper(e,F(1,3),F(1,2),F(1))
    terms['U4']=upper(a,mu,F(1,2),F(1,2))+upper(mu,d,F(57,64),F(5,2))
    def integrand(t,u,enhance):
        V=cap(t)
        if enhance:V=V.maximum(cap(t+u)).maximum(F(153,256)-F(3,8)*a)
        residual=V-t-u;s=residual/I(a)
        # Symbolic pointwise alpha>=2 proof; no dependency loss at equality.
        return 4*(sf.at(s,'b')+df.lookup(s,True))/(t*u*residual)
    step=(h-a)/pair_n;L5=I(0)
    physical=[I(a+i*step,a+(i+1)*step) for i in range(pair_n)]
    for i,t in enumerate(physical):
        for j in range(i):L5+=integrand(t,physical[j],True)*I(step*step)
        # The diagonal triangle has exact area step^2/2; the enclosing square
        # bounds its integrand, with the analytic insertion used only on u<t.
        L5+=integrand(t,t,True)*I(step*step/2)
    terms['L5']=L5
    print(json.dumps({'stage':'physical_L5','seconds':time.monotonic()-start}),flush=True)
    L6=I(0)
    for left,right,enhance in [(h,mu,True),(mu,d,False)]:
        dt=(right-left)/pair_n
        for i in range(pair_n):
            t=I(left+i*dt,left+(i+1)*dt)
            for u in physical:L6+=integrand(t,u,enhance)*I(dt*step)
    terms['L6']=L6
    print(json.dumps({'stage':'physical_L6','seconds':time.monotonic()-start}),flush=True)
    terms['U7']=K*boxes1(lambda t:(log((1-2*t)/t)+sf.at(1/t-1)-1)/(t*(1-t)),d,F(1,3),line_n)
    terms['U8']=K*boxes1(lambda t:log(2-3*t)/(t*(1-t)),a,F(1,3),line_n)
    terms['U9']=K*boxes1(lambda t:log((1-t-d)/I(d))/(t*(1-t)),h,d,line_n)
    terms['U10']=K*C/2*boxes1(lambda z:log(I(h)/z)**2*log(z/I(a))/z**2,a,h,line_n)
    terms['U11full']=K*C*log(I(d/h))*boxes1(lambda z:log(I(h)/z)*log(z/I(a))/z**2,a,h,line_n)
    c=(3*terms['L1']+terms['L2']+terms['L5']+terms['L6']-terms['U3']-terms['U4']-2*terms['U7']-terms['U8']-terms['U9']-terms['U10']-terms['U11full'])/4
    out={'status':'CERTIFIED_COMBINED_56_GENERATION_UNIFORM_SWITCH_FUNCTIONAL_UNDER_EXPLICIT_ANALYTIC_INPUTS',
         'asserted_lower':display(c),'safe_lower_1_90_verified':c.lo*10>19*Q,'safe_lower_1_907_verified':c.lo*1000>1907*Q,
         'parameters':dict(a=str(a),h=str(h),d=str(d),K=str(K)),'meshes':dict(pair=pair_n,line=line_n,ordinary=g),
         'terms':{n:display(x) for n,x in terms.items()},
         'inherited_audited_numerical_design':['Ordinary normalized A/a DDE rebuilt; shared primitive I/log only.',
             'L2 and U3/U4 correction integrals split exactly at rational table-cell inverse images.',
             'L5 covered by physical rectangles and diagonal half-squares; L6 by physical rectangles.',
             'Global symbolic alpha>=2 for a>=25/352 prevents cancellation losses on equality branches.',
             'All switched integrals and the eleven signed multiplicities rebuilt from their formulas.'],
         'limitations':['Lower endpoint alone bounds the actual functional from below; upper endpoint is not a true-functional upper bound.',
             'Analytic distribution, finite-table seed/closure and insertion are inherited separately, not established by this script.',
             'Shared I/log backend is explicitly not an independent low-level arithmetic implementation.',
             'No new theorem or imbalance exponent asserted by this numeric replay alone.'],
         'elapsed_seconds':time.monotonic()-start,
         'bindings_sha256':{n:sha(n) for n in [Path(__file__).name,'certify_conservative_chen.py','certify_combined_finite_table.json','certify_taylor_fifth_seed_transfer.json','certify_second_seed_omega_envelope.json','independent_finite_chen_functional.py','finite-chen-functional-independent-audit.json','combined-finite-table-independent-audit.json','finite-table-chen-insertion.tex','finite-w2-action-bank.tex']}}
    assert out['safe_lower_1_90_verified']
    (BASE/'combined-chen-functional-uniform.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'coefficient':display(c),'safe_1_907':out['safe_lower_1_907_verified'],'seconds':out['elapsed_seconds']}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pairs',type=int,default=384);ap.add_argument('--lines',type=int,default=8192);ap.add_argument('--ordinary',type=int,default=8192)
    x=ap.parse_args();run(x.pairs,x.lines,x.ordinary)
