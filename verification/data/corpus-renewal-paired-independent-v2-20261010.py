"""Independent enclosure of the same-mask paired refinement.

No project numerical implementation is imported.  Old floors are bounded
ABOVE using exact ordinary log endpoints, the literal chord geometry and
closed-cell h minima.  The additional positive mass uses a separate decimal
integer backend and alternating Mercator logarithms, not the pilot's atanh
backend.  This script never edits the manuscript or research state.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from math import ceil, floor
import hashlib
import json
import time

B = Path(__file__).resolve().parent
S = 10**80
TERMS = 320
STEM = 'corpus-renewal-paired-independent-v2-20261010'


class J:
    def __init__(self, x=0, hi=None, raw=False):
        if raw:
            self.lo, self.hi = int(x), int(hi)
        else:
            x = F(x)
            self.lo, self.hi = floor(x*S), ceil(x*S)
    def __add__(self, y):
        y = y if isinstance(y,J) else J(y)
        return J(self.lo+y.lo,self.hi+y.hi,raw=True)
    __radd__ = __add__
    def __neg__(self):
        return J(-self.hi,-self.lo,raw=True)
    def __sub__(self,y):
        return self+-asj(y)
    def __rsub__(self,y):
        return asj(y)+-self
    def __mul__(self,y):
        y=asj(y)
        a=[self.lo*y.lo,self.lo*y.hi,self.hi*y.lo,self.hi*y.hi]
        return J(min(a)//S,-(-max(a)//S),raw=True)
    __rmul__=__mul__
    def __truediv__(self,y):
        y=asj(y)
        assert y.lo>0
        a=[(self.lo*S,y.lo),(self.lo*S,y.hi),(self.hi*S,y.lo),(self.hi*S,y.hi)]
        return J(min(n//d for n,d in a),max(-(-n//d) for n,d in a),raw=True)
    def __rtruediv__(self,y):
        return asj(y)/self


def asj(y):
    return y if isinstance(y,J) else J(y)


@lru_cache(None)
def mercator(q):
    q=F(q)
    assert 1<=q<=F(3,2)
    y=J(q-1); p=J(1); s=J(0)
    # Even partial sum is below log(1+y); next term bounds its error.
    for k in range(1,TERMS+1):
        p=p*y
        s=s+(p/k if k%2 else -(p/k))
    tail=(p*y)/(TERMS+1)
    assert tail.lo>=0
    return J(s.lo,s.hi+tail.hi,raw=True)


@lru_cache(None)
def logq(q):
    q=F(q)
    assert 1<=q<=2
    shifts=0
    while q>F(3,2):
        q*=F(4,5);shifts+=1
    return mercator(q)+shifts*mercator(F(5,4))


def pack(x):
    return dict(lower=str(F(x.lo,S)),upper=str(F(x.hi,S)),
                display_lower=x.lo/S,display_upper=x.hi/S)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    start=time.time()
    names=['positive-pair-deletion-certificate-20261010.json',
           'paired-y6-y11-zero-strip-certificate-20261010.json',
           'expanded-w2-midpoint-table248-v1.json',
           'certify_sharp_184_baseline.py','certify_combined_chen_functional.py',
           'certify_extended_prefix.py','sharp-248-pair-refinement-512.json',
           'positive-pair-deletion-selected-parent-cover-20261010.json',
           'paired-y6-y11-parent-cover-root-check-20261010.json',
           'corpus-renewal-omega-independent-20261010.txt',
           'corpus-renewal-omega-independent-20261010.json']
    bindings={n:sha(B/n) for n in names}
    pos=json.loads((B/names[0]).read_text())
    zero=json.loads((B/names[1]).read_text())
    bank=json.loads((B/names[2]).read_text())
    e1=[tuple(map(int,p)) for p in pos['selected_pairs']]
    e0=[(r['t_index'],j) for r in zero['strips'] for j in range(r['first_u_index'],512)]
    assert len(e1)==len(set(e1))==11481 and len(e0)==len(set(e0))==6091
    assert not set(e0)&set(e1)
    mask=''.join(f'{i},{j}\n' for i,j in e1).encode()
    maskhash=hashlib.sha256(mask).hexdigest()
    assert maskhash==pos['selected_mask_sha256']
    a,h,mu,d,C,K=map(F,['417/5000','19/125','25/128','59/250','57/64','8452171/1250000'])
    dt,du=(d-mu)/512,(h-a)/512
    hh=list(map(int,bank['h_cell_lower_numerators']))
    Q=int(bank['fixed_point_denominator']); grid=int(bank['grid_denominator'])
    assert Q==2**80 and grid==50
    extras={};lost=J(0);rows=[];cross=0;minimum_inner=F(1)
    for group,cells in [('E0',e0),('E1',e1)]:
        extra=J(0)
        for i,j in cells:
            tl,tr=mu+i*dt,mu+(i+1)*dt
            ul,ur=a+j*du,a+(j+1)*du
            lt,lu=logq(tr/tl),logq(ur/ul);mass=lt*lu
            zl=min(ul,F(3,13)*(1-tr-ur))
            assert a<zl<=ul and 1-tr-ur-2*zl>=F(7,3)*zl
            minimum_inner=min(minimum_inner,zl)
            fibre=1/J(a)-(1+logq(zl/a))/J(zl)
            assert fibre.lo>0
            gross=K/F(80)*mass*fibre
            extra+=gross
            rmin,rmax=C-F(5,2)*tr-ur,C-F(5,2)*tl-ul
            sl,sh=rmin/a,rmax/a
            if group=='E0':
                assert sh<=2
                continue
            assert 2<sh<=F(9,4)
            if sl<2:
                upper=J(0);cross+=1
            else:
                fl,fh=logq(sl-1)/J(sl),logq(sh-1)/J(sh)
                ordinary_constant=fl*mass/J(a)
                moment=C*mass-F(5,2)*(tr-tl)*lu-(ur-ul)*lt
                left=rmax*mass-moment;right=moment-rmin*mass
                assert left.lo>0 and right.lo>0
                chord=(fl*left+fh*right)/J(a*(rmax-rmin))
                ordinary=J(max(ordinary_constant.lo,chord.lo),max(ordinary_constant.hi,chord.hi),raw=True)
                # All closed table cells meeting [sl,sh].  The actual outward
                # lookup can only meet MORE cells, so its minimum is no larger.
                first=max(0,ceil((sl-1)*grid)-1)
                last=min(len(hh)-1,floor((sh-1)*grid))
                hc=F(min(hh[first:last+1]),Q)
                upper=ordinary+hc*mass/J(rmax)
            lost+=upper
            rows.append([i,j,str(upper.hi),str(gross.lo),str(gross.hi)])
        extras[group]=pack(extra)
    assert cross==227
    gross=J(F(extras['E0']['lower']))+J(F(extras['E1']['lower']))
    # Above gross is a rigorous lower enclosure of the sum of the two lowers.
    safe_gross=F(558583,10**9)
    safe_loss=F(466977,200000000)
    old_loss=F(1186047,500000000)
    assert F(gross.lo,S)>safe_gross
    assert F(lost.hi,S)<safe_loss
    assert safe_loss<old_loss
    safe_gain=safe_gross+old_loss-safe_loss
    assert safe_gain==F(595792,10**9)
    for n,digest in bindings.items():assert sha(B/n)==digest,n
    bindings[Path(__file__).name]=sha(Path(__file__))
    out=dict(status='INDEPENDENT_DIRECTED_ENCLOSURES_ON_IDENTICAL_MASK',
      date='2026-10-10',buchstab_threshold='7/3',total_cells=len(e0)+len(e1),backend='decimal integer intervals; alternating Mercator; no project backend import',
      scale=str(S),log_terms=TERMS,selected_mask_sha256=maskhash,
      cell_counts=dict(E0=len(e0),E1=len(e1),E1_zero_crossing=cross),
      extra_gross=extras,discarded_author_floor_upper_proxy=pack(lost),
      safe_additional_gross=str(safe_gross),safe_E1_loss_upper=str(safe_loss),
      old_E1_loss=str(old_loss),safe_loss_recovery=str(old_loss-safe_loss),
      safe_raw_increment=str(safe_gain),additional_future_slope_deduction=str(safe_gross/K),
      minimum_inner_Z=str(minimum_inner),
      loss_semantics='Upper bound on exact saved author floors, not on actual arithmetic Y6 count. Log endpoints dominate certified original sf lower endpoints.',
      row_columns=['t_index','u_index','loss_upper_S','extra_gross_lower_S','extra_gross_upper_S'],
      rows=rows,bindings_sha256=bindings,seconds=time.time()-start)
    dest=B/(STEM+'.json')
    assert not dest.exists(),'immutable independent result already exists'
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:out[k] for k in ['status','extra_gross','discarded_author_floor_upper_proxy','safe_raw_increment','seconds']},indent=2))


if __name__=='__main__':main()
