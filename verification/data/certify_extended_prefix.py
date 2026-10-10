"""Exact enclosure of the extended first-stop, variable third-prefix functional.

NUMERICAL ONLY: this program does not prove the analytic sieve interface.
Every bound and subdivision priority uses integers or Fractions. Floats are
used only in progress and final display fields.

The normalized linear-sieve functions are fbar=f/(2 exp(gamma)) and
Fbar=F/(2 exp(gamma)). Standard positivity and monotonicity of f,F are inputs.
L=s*fbar and H=s*Fbar obey L'=Fbar(s-1), H'=fbar(s-1).
For s>=1, Fbar'' is nonnegative and <=2, piecewise at the integer breaks;
for s>=2, fbar'' is nonpositive with absolute value <=1. Indeed Fbar<=1,
|Fbar'|<=1/s and 0<=fbar'<=1/s; differentiating their delay equations gives
these looser global bounds. Thus trapezoidal integration over a grid cell
gives a lower/upper enclosure with error <= h^3/6 for L and <= h^3/12 for H.
All integer breaks are grid points. L is concave with |L''|<=1 for s>=2,
so linear interpolation has upper error <=h^2/8. H is convex with H''<=1.

For D=Fbar-fbar, D decreases and (sD)'=-D(s-1), s>=3. Hence
D(s)<=D(s0)*(s0/s)^2 for s>=s0>=3. Also fbar'(s)<=D(s-1)/s for s>=3.
Both give additional safe bounds for differences; no near-equal floats
are subtracted. All above are mathematical inputs/proofs, not empirical fits.
"""
import argparse
import hashlib
import heapq
import json
import time
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path

from certify_conservative_chen import Q, BITS, I, log, ceildiv

THETA = F(299, 500)
MU = F(25, 128)
RHO = F(19, 96)
TAIL_CUTOFF = F(1, 100)
TAIL_S0 = 10
DEN = 48000 * 2**24


def coord(q):
    v = q*DEN
    assert v.denominator == 1
    return v.numerator


NT, NMU, NRHO, NK = map(coord, (THETA, MU, RHO, TAIL_CUTOFF))


class SieveTable:
    """Fixed-point trapezoidal enclosures, with rigorous curvature errors."""
    def __init__(self, n=16384, maximum=12):
        assert n > 0 and maximum >= 4
        self.n, self.maximum = n, maximum
        size = maximum*n+1
        self.ll, self.lh = [0]*size, [0]*size
        self.hl, self.hh = [Q]*size, [Q]*size
        err_l = ceildiv(Q, 6*n**3)
        err_h = ceildiv(Q, 12*n**3)
        for i in range(2*n+1, size):
            j = i-n
            # L'=Fbar(s-1), which is convex. Trapezoid is an upper bound.
            lo = (self.hl[j-1]*n//(j-1) + self.hl[j]*n//j)//(2*n)
            hi = ceildiv(ceildiv(self.hh[j-1]*n,j-1)
                         + ceildiv(self.hh[j]*n,j), 2*n)
            self.ll[i] = max(self.ll[i-1], self.ll[i-1]+lo-err_l)
            self.lh[i] = self.lh[i-1]+hi
            if i > 3*n:
                # H'=fbar(s-1), concave on this range. Trapezoid is lower.
                lo = (self.ll[j-1]*n//(j-1) + self.ll[j]*n//j)//(2*n)
                hi = ceildiv(ceildiv(self.lh[j-1]*n,j-1)
                             + ceildiv(self.lh[j]*n,j), 2*n)
                self.hl[i] = self.hl[i-1]+lo
                self.hh[i] = self.hh[i-1]+hi+err_h
        self.interp_error = ceildiv(Q,8*n*n)
        self.d_integer = {
            s: ceildiv(self.hh[s*n]-self.ll[s*n],s)
            for s in range(3,maximum+1)
        }
        assert all(v>0 for v in self.d_integer.values())
        for s in (F(5,2),F(3),F(7,2),F(4)):
            lo,hi=self.point(s.numerator,s.denominator)
            exact=log(I(s-1))/I(s)
            assert lo<=exact.lo<=exact.hi<=hi

    def numerator(self, a, b, lower, upper, convex=False):
        """Enclose L(a/b) or H(a/b) by interpolation plus curvature."""
        q,r=divmod(a*self.n,b)
        assert q<=self.maximum*self.n
        if r==0:
            return lower[q],upper[q]
        lo=(lower[q]*(b-r)+lower[q+1]*r)//b
        hi=ceildiv(upper[q]*(b-r)+upper[q+1]*r,b)
        if convex:
            lo-=self.interp_error
        else:
            hi+=self.interp_error
        return lo,hi

    def point(self, a, b):
        if a<=2*b:
            return 0,0
        if a>self.maximum*b:
            # f increases; F decreases and both converge to a common limit.
            m=self.maximum
            return self.ll[m*self.n]//m,ceildiv(self.hh[m*self.n],m)
        lo,hi=self.numerator(a,b,self.ll,self.lh)
        return lo*b//a,ceildiv(hi*b,a)

    def d_upper(self,a,b):
        """Upper bound for D(a/b), only called when a/b>=3."""
        assert a>=3*b
        last=min(self.maximum,a//b)
        return min(ceildiv(self.d_integer[s]*s*s*b*b,a*a)
                   for s in range(3,last+1))


def v_lo(lo,hi):
    return min((DEN+lo)//2,(57*DEN-96*hi)//64)


def v_point_hi(t):
    return min(ceildiv(DEN+t,2),ceildiv(57*DEN-96*t,64))


def v_hi(lo,hi):
    if lo<=NMU<=hi:
        return ceildiv(153*DEN,256)
    return max(v_point_hi(lo),v_point_hi(hi))


def phi_range(cell):
    """Range enclosure of the exact reduced full envelope on ordered points.

    For 0<=z<=y<=x<=mu, the middle affine piece of v is redundant,
    w1's inner minimum is B(z)=153/256-3*z/8, and w2 is dominated.
    Each pair gate is automatically satisfied by at least 3/256 whenever
    v(pair)>=v(x): for pair s<=mu use 1/2-5*mu/2; for s>=mu,
    v(s)>=v(x) implies s<=(25/32-x)/3 and hence
    v(s)-(s+2*x)>=23/96-7*x/6>=3/256. Thus pair gates can be removed
    inside the max with v(x). Formula bounds are evaluated on a rectangle,
    but are needed only on its ordered subset.
    """
    x0,x1,y0,y1,z0,z1=cell
    # Initial x intervals are split at mu and subdivision preserves this.
    # On an ordered point with x>mu, use only v(x). On a symmetric square
    # above mu, extend it as v(max(x,y)); its range is the same interval.
    assert x1<=NMU or x0>=NMU
    if x0>=NMU:
        return min(NT,v_lo(x0,x1)),min(NT,v_hi(x0,x1))
    branches=[(x0,x1),(y0,y1),(x0+y0,x1+y1),(x0+y0+z0,x1+y1+z1),
              (x0+z0,x1+z1),(y0+z0,y1+z1)]
    low=max(v_lo(a,b) for a,b in branches)
    high=max(v_hi(a,b) for a,b in branches)
    low=max(low,(153*DEN-96*z1)//256)
    high=max(high,ceildiv(153*DEN-96*z0,256))
    return min(NT,low),min(NT,high)


@lru_cache(maxsize=500000)
def logarithmic_width(lo,hi):
    ans=log(I(F(hi,lo)))
    return ans.lo,ans.hi


def density_mass(cell):
    """Integral of 4/(x*y*z^2) over the complete rectangle, in Q units."""
    x0,x1,y0,y1,z0,z1=cell
    a,b=logarithmic_width(x0,x1),logarithmic_width(y0,y1)
    lo=a[0]*b[0]//Q
    hi=ceildiv(a[1]*b[1],Q)
    numerator=4*DEN*(z1-z0)
    denominator=z0*z1
    return lo*numerator//denominator,ceildiv(hi*numerator,denominator)


def affine_extreme(cell,p,maximum):
    a,b,c,d,e=p
    endpoints=cell[1::2] if maximum else cell[::2]
    opposite=cell[::2] if maximum else cell[1::2]
    return a*DEN+sum(co*(hi if co>=0 else lo)
                     for co,hi,lo in zip((b,c,d),endpoints,opposite))


def plane_difference(p,q):
    return tuple(p[i]*q[4]-q[i]*p[4] for i in range(4))+(1,)


def active_min_plane(cell,planes):
    for p in planes:
        if all(affine_extreme(cell,plane_difference(p,q),True)<=0 for q in planes):
            return p
    return None


def v_planes(b,c,d):
    return [(1,b,c,d,2),(57,-96*b,-96*c,-96*d,64)]


def affine_cap_lower(cell):
    """Return one whole-box affine minorant of min(Theta,Phi), if available."""
    if cell[0]>=NMU:
        # Used only for wholly ordered rectangles; a symmetric square above
        # mu requires the separate exact mean of max(x,y) below.
        return (57,-96,0,0,64)
    choices=[]
    for b,c,d in ((1,0,0),(0,1,0),(1,1,0),(1,1,1),(1,0,1),(0,1,1)):
        p=active_min_plane(cell,v_planes(b,c,d))
        if p is not None:
            choices.append(p)
    choices.append((153,0,0,-96,256))
    choices=[p for p in choices if affine_extreme(cell,p,True)<=NT*p[4]]
    if not choices:
        return None
    midpoint=tuple((cell[2*i]+cell[2*i+1])//2 for i in range(3))
    def at_mid(p):
        return p[0]*DEN+sum(a*b for a,b in zip(p[1:4],midpoint))
    best=choices[0]
    for p in choices[1:]:
        if at_mid(p)*best[4]>at_mid(best)*p[4]:
            best=p
    return best


def jensen_difference_upper(cell,table):
    """Weighted integral mean bound, valid only on a full raw-domain cell.

    The density factors across x,y,z. Jensen bounds the first concave fbar
    from above. An endpoint chord bounds the second concave fbar from below.
    The selected affine cap is a whole-cell minorant of the actual cap.
    This removes first-order argument variation on cells inside affine pieces.
    """
    x0,x1,y0,y1,z0,z1=cell
    above_symmetric=(x0>=NMU and x0==y0 and x1==y1 and z1<=y0)
    p=affine_cap_lower(cell)
    if p is None:
        return None
    a,b,c,d,e=p
    # Every candidate has b/e,c/e <=1 and a positive residual numerator.
    if b>e or c>e or a*DEN+(b-e)*x1+(c-e)*y1<=0:
        return None
    old_lo=F(a*DEN+(b-e)*x1+(c-e)*y1+(d-e)*z1,e*z1)
    old_hi=F(a*DEN+(b-e)*x0+(c-e)*y0+(d-e)*z0,e*z0)
    if old_lo<2<old_hi:
        return None
    lx=I(*logarithmic_width(x0,x1),raw=True)
    ly=I(*logarithmic_width(y0,y1),raw=True)
    ex=I(x1-x0)/lx
    ey=I(y1-y0)/ly
    ez=I(F(z0+z1,2*z0*z1))
    new_mean=(I(NT)-ex-ey)*ez-1
    new_hi=table.point(new_mean.hi,Q)[1]
    if old_hi<=2:
        return new_hi
    if above_symmetric:
        # Under independent log-uniform x,y on [a,b],
        # E max(x,y)=2*(b*log(b/a)-b+a)/log(b/a)^2.
        # Therefore this symmetric extension has the exact old mean below.
        # Endpoint extremes remain those at (x0,y0,z0),(x1,y1,z1).
        emax=2*(x1*lx-x1+x0)/(lx*lx)
        old_mean=(I(F(57*DEN,64))-ex-ey-F(3,2)*emax)*ez-1
    else:
        old_mean=(I(F(a*DEN,e))+F(b-e,e)*ex+F(c-e,e)*ey)*ez+F(d-e,e)
    fl=table.point(old_lo.numerator,old_lo.denominator)[0]
    fh=table.point(old_hi.numerator,old_hi.denominator)[0]
    if old_lo==old_hi:
        old_lower=fl
    else:
        weight=(old_mean-I(old_lo))/I(old_hi-old_lo)
        weight=I(max(0,weight.lo),min(Q,weight.hi),raw=True)
        chord=I(fl,fl,raw=True)+weight*I(fh-fl,fh-fl,raw=True)
        old_lower=chord.lo
    return max(0,new_hi-old_lower)


def enclosure(cell, table):
    x0,x1,y0,y1,z0,z1=cell
    if y0>x1 or z0>y1 or x0+y0+3*z0>NT:
        return 0,0,'outside'
    # For max(x,y)<=mu use the symmetric reduced Phi. Above mu extend the
    # ordered cap by v(max(x,y)), which is also symmetric. Initial x/y cuts
    # at mu ensure an identical-interval square never crosses the jump.
    # In either case the ordered half has half the full-square integral.
    symmetric_half=(x0==y0 and x1==y1 and z1<=y0)
    inside=((y1<=x0 or symmetric_half) and z1<=y0 and x1+y1+3*z1<=NT)
    cap_lo,cap_hi=phi_range(cell)
    if cap_lo>=NT:
        return 0,0,'full_cap'
    # Both levels minus x-y are positive throughout this compact cube.
    assert cap_lo-x1-y1>0
    nlo=NT-x1-y1-z1
    nhi=NT-x0-y0-z0
    olo=cap_lo-x1-y1-z1
    ohi=cap_hi-x0-y0-z0
    fnlo=table.point(nlo,z1)[0]
    fnhi=table.point(nhi,z0)[1]
    folo=table.point(olo,z1)[0]
    fohi=table.point(ohi,z0)[1]
    dl=max(0,fnlo-fohi) if inside else 0
    dh=max(0,fnhi-folo)
    # fbar is increasing and concave above 2. For fixed delta>=0,
    # fbar(s+delta)-fbar(s) decreases when s>=2, and its global maximum
    # over s<=2 occurs at s=2. This keeps the common argument variation
    # from being counted twice in a difference of near-equal functions.
    sa=max(2*z1,olo)
    delta=NT-cap_lo
    paired_hi=table.point(sa*z0+delta*z1,z1*z0)[1]
    paired_lo=table.point(sa,z1)[0]
    dh=min(dh,max(0,paired_hi-paired_lo))
    if inside and olo>=2*z1:
        # The same monotonic difference argument gives a correlated lower
        # bound using largest possible old argument and smallest increment.
        paired_low=table.point(ohi*z1+(NT-cap_hi)*z0,z0*z1)[0]
        paired_high=table.point(ohi,z0)[1]
        dl=max(dl,max(0,paired_low-paired_high))
    if olo>=3*z1:
        dh=min(dh,table.d_upper(olo,z1))
    if olo>=4*z1:
        # f'(s)<=D(s-1)/s and s_new-s_old <=(Theta-cap_lo)/z0.
        dh=min(dh,ceildiv(table.d_upper(olo-z1,z1)*(NT-cap_lo)*z1,
                          z0*olo))
    if dh==0:
        return 0,0,'zero_loss'
    if inside:
        jensen=jensen_difference_upper(cell,table)
        if jensen is not None:
            dh=min(dh,jensen)
    ml,mh=density_mass(cell)
    if symmetric_half:
        ml//=2
        mh=ceildiv(mh,2)
    lower=dl*ml//Q
    upper=ceildiv(dh*mh,Q)
    assert 0<=lower<=upper,(cell,lower,upper)
    return lower,upper,'interior' if inside else 'boundary'


def first_integral(terms):
    """(4/Theta)*integral_2^(Theta/rho-1) log(s-1)/s ds."""
    # q=s-2. log(1+q)/(2+q) has alternating coefficients a_k>0,
    # a_1=1/2, a_k=a_(k-1)/2+1/(2*k). Induction gives
    # a_k>=1/(k+1) and hence a_(k+1)<=a_k. Since q_max<1,
    # the integrated terms decrease: even partial sums are lower,
    # odd partial sums upper. Use consecutive even/odd truncations.
    assert terms>=2 and terms%2==0
    q=THETA/RHO-3
    assert 0<q<1
    coefficient=F(1,2);partial=F(0)
    for k in range(1,terms+2):
        if k>1:
            coefficient=coefficient/2+F(1,2*k)
        term=coefficient*q**(k+1)/(k+1)
        partial+=term if k%2 else -term
        if k==terms:
            low=4*partial/THETA
    high=4*partial/THETA
    return low.numerator*Q//low.denominator,ceildiv(high.numerator*Q,high.denominator)


def initial_cells():
    # Geometric t3 strips, and exact x/y cuts at the discontinuity mu.
    cuts=[NK]
    while cuts[-1]<NRHO:
        cuts.append(min(NRHO,2*cuts[-1]))
    cuts=sorted(set(cuts+[NMU]))
    for z0,z1 in zip(cuts,cuts[1:]):
        axis_cuts=[z0]+([NMU] if z0<NMU else [])+[NRHO]
        for x0,x1 in zip(axis_cuts,axis_cuts[1:]):
            for y0,y1 in zip(axis_cuts,axis_cuts[1:]):
                yield (x0,x1,y0,y1,z0,z1)


def tail_upper(table):
    """Upper bound for the entire ordered 0<z<h third correction.

    Every affine branch of v(x)-2x is >=rho for 0<=x<=rho, while
    Theta-2x>=rho. Thus L-S>=rho-z. For h<=rho/(s0+1),
    s_old>=s0 and loss<=D(s_old)<=d0*s0^2*z^2/(rho-h)^2.
    Integrating 4/(x*y*z^2) on z<=y<=x<=rho leaves
    2*d0*s0^2/(rho-h)^2 * integral_0^h log(rho/z)^2 dz.
    All operations below are outward integer interval operations.
    """
    assert 3<=TAIL_S0<=table.maximum
    assert 0<TAIL_CUTOFF<=RHO/(TAIL_S0+1)
    d0=I(table.d_integer[TAIL_S0],table.d_integer[TAIL_S0],raw=True)
    ell=log(I(RHO/TAIL_CUTOFF))
    upper=2*d0*TAIL_S0**2*TAIL_CUTOFF*(ell*ell+2*ell+2)/I((RHO-TAIL_CUTOFF)**2)
    return upper.hi


def main(args):
    assert MU<RHO and 3*RHO<THETA<F(3,5)
    # The fallback v(x) above mu is everywhere below this global cap.
    assert THETA>F(153,256)
    start=time.monotonic()
    table=SieveTable(args.table_mesh,args.table_maximum)
    first_lo,first_hi=first_integral(args.first_terms)
    tail_hi=tail_upper(table)
    fixed_lo=4/THETA+F(first_lo,Q)
    fixed_hi=4/THETA+F(first_hi+tail_hi,Q)
    heap=[];serial=0;total_lo=total_hi=0
    counts={}
    for cell in initial_cells():
        lo,hi,kind=enclosure(cell,table)
        total_lo+=lo;total_hi+=hi
        if hi>lo:
            serial+=1
            heapq.heappush(heap,(-(hi-lo),serial,cell,lo,hi))
    target=F(args.target)
    splits=0
    def progress():
        return {'splits':splits,'cells':len(heap),
                'K_lower_display_only':float(fixed_lo+F(total_lo,Q)),
                'K_upper_display_only':float(fixed_hi+F(total_hi,Q)),
                'seconds':round(time.monotonic()-start,2)}
    print(json.dumps(progress()),flush=True)
    while splits<args.max_splits and fixed_hi+F(total_hi,Q)>=target:
        _,_,cell,lo,hi=heapq.heappop(heap)
        total_lo-=lo;total_hi-=hi
        # Relative widths improve resolution near the small third cutoff.
        axis=0
        for candidate in (1,2):
            if ((cell[2*candidate+1]-cell[2*candidate])*cell[2*axis]
                    >(cell[2*axis+1]-cell[2*axis])*cell[2*candidate]):
                axis=candidate
        j=2*axis;mid=(cell[j]+cell[j+1])//2
        assert cell[j]<mid<cell[j+1]
        for side in (0,1):
            child=list(cell);child[j+1 if side==0 else j]=mid;child=tuple(child)
            cl,ch,kind=enclosure(child,table)
            counts[kind]=counts.get(kind,0)+1
            total_lo+=cl;total_hi+=ch
            if ch>cl:
                serial+=1
                heapq.heappush(heap,(-(ch-cl),serial,child,cl,ch))
        splits+=1
        if splits%args.report_every==0:
            print(json.dumps(progress()),flush=True)
    upper=fixed_hi+F(total_hi,Q)
    lower=fixed_lo+F(total_lo,Q)
    here=Path(__file__)
    dependency=here.with_name('certify_conservative_chen.py')
    result={'status':'NUMERICAL_ONLY_EXTENDED_PREFIX_ANALYTIC_INTERFACE_SEPARATE',
            'theta':str(THETA),'first_stop':str(RHO),'three_prime_mu':str(MU),
            'compact_t3_lower':str(TAIL_CUTOFF),'tail_s0':TAIL_S0,
            'target_K':str(target),'strict_upper_target_verified':upper<target,
            'lower_exact':str(lower),'upper_exact':str(upper),
            'lower_display_only':float(lower),'upper_display_only':float(upper),
            'gap_below_target':str(target-upper),
            'first_integral_lower_exact':str(F(first_lo,Q)),
            'first_integral_upper_exact':str(F(first_hi,Q)),
            'third_tail_lower_exact':'0','third_tail_upper_exact':str(F(tail_hi,Q)),
            'third_integral_lower_exact':str(F(total_lo,Q)),
            'third_integral_upper_exact':str(F(total_hi,Q)),
            'fixed_point_bits':BITS,'coordinate_denominator':DEN,
            'sieve_table_mesh':args.table_mesh,'sieve_table_maximum':args.table_maximum,
            'first_integral_even_series_terms':args.first_terms,'splits':splits,
            'remaining_cells':len(heap),'counts':counts,
            'D_integer_upper_exact':{str(s):str(F(v,Q)) for s,v in table.d_integer.items()},
            'seconds':time.monotonic()-start,
            'source_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
            'arithmetic_dependency_sha256':hashlib.sha256(dependency.read_bytes()).hexdigest(),
            'scope':['Encloses only the displayed first-stop plus variable-third-cap integral.',
                     'Uses standard positivity and monotonicity of the normalized linear-sieve functions.',
                     'Uses reduced full Phi for t1<=25/128 and only v(t1) above it up to 19/96.',
                     'The noncompact region below t3=1/100 has a certified D(s) tail upper bound.',
                     'The extension across x/y order uses Phi below mu and v(max(x,y)) above mu.',
                     'Does not prove distribution, factorization, raw-density asymptotics or a counting theorem.']}
    out=Path(args.output) if args.output else here.with_suffix('.json')
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--target',default='681/100')
    ap.add_argument('--max-splits',type=int,default=200000)
    ap.add_argument('--table-mesh',type=int,default=16384)
    ap.add_argument('--table-maximum',type=int,default=12)
    ap.add_argument('--first-terms',type=int,default=20)
    ap.add_argument('--report-every',type=int,default=10000)
    ap.add_argument('--output',default=None)
    main(ap.parse_args())
