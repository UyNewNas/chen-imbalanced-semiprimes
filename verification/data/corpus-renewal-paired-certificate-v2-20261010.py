"""Directed rational certificate for refining an already adopted paired mask.

Only additional gross mass and the replacement saved-floor loss are used.
The same old E0/E1 mask, coefficients, table and parent ownership survive.
No old output, state, main, formalization or selected parameter is changed.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib, json, time

B=Path(__file__).resolve().parent
P=1<<160; W=1<<224
a,h,d,mu,C,K=F(417,5000),F(19,125),F(59,250),F(25,128),F(57,64),F(8452171,1250000)
def floor(q):q=F(q);return q.numerator//q.denominator
def ceil(q):q=F(q);return -((-q.numerator)//q.denominator)
def read(n):return json.loads((B/n).read_text(encoding='utf-8-sig'))
def path(n):
    if n.startswith('/') and len(n)>2 and n[2]==':':n=n[1:]
    p=Path(n);return p if p.is_absolute() else B/p
def sha(n):
    with path(n).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

class Interval:
    def __init__(self,x=0,y=None,raw=False):
        if raw:self.lo,self.hi=int(x),int(y)
        else:self.lo,self.hi=floor(F(x)*P),ceil(F(x)*P)
        assert self.lo<=self.hi
    def __add__(self,q):
        if not isinstance(q,Interval):q=Interval(q)
        return Interval(self.lo+q.lo,self.hi+q.hi,True)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo,True)
    def __sub__(self,q):
        if not isinstance(q,Interval):q=Interval(q)
        return self+-q
    def __mul__(self,q):
        if not isinstance(q,Interval):q=Interval(q)
        products=[self.lo*q.lo,self.lo*q.hi,self.hi*q.lo,self.hi*q.hi]
        return Interval(min(products)//P,ceil(F(max(products),P)),True)
    __rmul__=__mul__
    def __truediv__(self,q):
        if not isinstance(q,Interval):q=Interval(q)
        assert q.lo>0
        lows=[n*P//v for n in (self.lo,self.hi) for v in (q.lo,q.hi)]
        highs=[ceil(F(n*P,v)) for n in (self.lo,self.hi) for v in (q.lo,q.hi)]
        return Interval(min(lows),max(highs),True)
    def record(self):return {'lower':str(F(self.lo,P)),'upper':str(F(self.hi,P))}

TAIL=ceil(F(9,4*193*3**193)*W)
@lru_cache(maxsize=None)
def log_interval(q):
    q=F(q);assert 1<=q<=2
    z=(q-1)/(q+1);zl=floor(z*W);zh=ceil(z*W)
    z2l=zl*zl//W;z2h=(zh*zh+W-1)//W
    tl,th=zl,zh;sl=sh=0
    for j in range(96):
        k=2*j+1;sl+=tl//k;sh+=(th+k-1)//k
        tl=tl*z2l//W;th=(th*z2h+W-1)//W
    shift=1<<64
    return Interval((2*sl)//shift,(2*sh+TAIL+shift-1)//shift,True)

BINDING_KEYS={'bindings_sha256','dependencies_sha256','verified_input_hashes','sources_sha256','sources_read_sha256'}
def named_bindings(obj):
    if isinstance(obj,dict):
        for key,value in obj.items():
            if key in BINDING_KEYS and isinstance(value,dict):
                for name,digest in value.items():
                    if isinstance(digest,str) and len(digest)==64:
                        yield name,digest
            elif isinstance(value,(list,dict)):
                yield from named_bindings(value)
    elif isinstance(obj,list):
        for value in obj:yield from named_bindings(value)

def mask_digest(mask):return hashlib.sha256(''.join(f'{i},{j}\n' for i,j in mask).encode()).hexdigest()

def main():
    start=time.monotonic();dest=B/'corpus-renewal-paired-certificate-v2-20261010.json'
    assert not dest.exists()
    direct=[Path(__file__).name,'corpus-renewal-paired-certificate-v2-20261010.txt',
      'positive-pair-deletion-certificate-20261010.json',
      'paired-y6-y11-zero-strip-certificate-20261010.json',
      'positive-pair-deletion-selected-parent-cover-20261010.json',
      'paired-y6-y11-parent-cover-root-check-20261010.json',
      'positive-pair-deletion-improvement-assembly-20261010.json',
      'positive-pair-deletion-improvement-20261010.txt',
      'paired-positive-cell-loss-interface-independent-20261010.txt',
      'chen-transfer-paired-deletion-independent-20261010.txt',
      'global-chen-lower-completion-20261010.txt',
      'corpus-renewal-paired-certificate-20261010.json',
      'sharp-248-pair-refinement-512.json',
      'independent-sharp-248-pair512-replay-20261008.json',
      'certify_sharp_248_pair_refinement.py','certify_sharp_184_baseline.py',
      'certify_conservative_chen.py','certify_extended_prefix.py',
      'certify_combined_chen_functional.py','expanded-w2-midpoint-table248-v1.json']
    bindings={n:sha(n) for n in direct}
    # Include and verify every directly named source/witness dependency in
    # the certificates supplying masks, floors, table, and ownership.
    roots=['positive-pair-deletion-certificate-20261010.json',
      'paired-y6-y11-zero-strip-certificate-20261010.json',
      'positive-pair-deletion-selected-parent-cover-20261010.json',
      'paired-y6-y11-parent-cover-root-check-20261010.json',
      'sharp-248-pair-refinement-512.json',
      'independent-sharp-248-pair512-replay-20261008.json',
      'expanded-w2-midpoint-table248-v1.json']
    source_maps_checked=0
    for name in roots:
        for child,digest in named_bindings(read(name)):
            assert child!='research-state.json'
            assert sha(child)==digest,(name,child)
            assert child not in bindings or bindings[child]==digest
            bindings[child]=digest;source_maps_checked+=1
    state_guard=sha('research-state.json');main_guard=sha('imbalanced-chen-transfer.tex')
    old=read('positive-pair-deletion-certificate-20261010.json')
    zold=read('paired-y6-y11-zero-strip-certificate-20261010.json')
    assembly=read('positive-pair-deletion-improvement-assembly-20261010.json')['selected_result']
    e1=[tuple(x) for x in old['selected_pairs']]
    e0=sorted((s['t_index'],j) for s in zold['strips'] for j in range(s['first_u_index'],512))
    assert len(e0)==6091 and len(e1)==11481 and len(set(e0+e1))==17572
    assert e0==sorted(set(e0)) and e1==sorted(set(e1))
    assert mask_digest(e1)==old['selected_mask_sha256']
    cover1=read('positive-pair-deletion-selected-parent-cover-20261010.json')
    cover0=read('paired-y6-y11-parent-cover-root-check-20261010.json')
    assert cover1['selected_mask_sha256']==mask_digest(e1)
    e1cover=sorted((s['t_index'],j) for s in cover1['complete_selected_pair_mask']
        for l,r in s['u_index_intervals_half_open'] for j in range(l,r))
    assert e1cover==e1
    assert cover1['methods']=={'0':61956} and cover1['priority40000_overlap']==0
    assert cover0['all_saved_methods_zero'] and cover0['priority40000_overlap']==0
    assert cover0['intersecting_parents']==34405
    assert F(assembly['safe_K'])==K
    # Literal deterministic old-cell replay. No new lower sieve estimate
    # is substituted for the floor whose deletion the old proof selected.
    import certify_sharp_184_baseline as author
    sf=author.Ordinary(16384);df=author.Defects('expanded-w2-midpoint-table248-v1.json');cache={}
    assert author.Q==1<<80
    dt=(d-mu)/512;du=(h-a)/512
    ts=[mu+i*dt for i in range(513)];us=[a+j*du for j in range(513)]
    lt=[log_interval(ts[i+1]/ts[i]) for i in range(512)]
    lu=[log_interval(us[j+1]/us[j]) for j in range(512)]
    la=[log_interval(x/a) for x in us]
    assert F(27,49)>F(11,20) and F(5,9)>F(11,20)
    assert log_interval(F(2)).lo*3>2*P
    assert F(1,2)+F(1,6)/F(10,3)==F(11,20)
    results={};rows=[];all_extra=Interval(0);lost=0;crossing=0
    for label,mask in [('E0',e0),('E1',e1)]:
        extra=Interval(0);base=Interval(0);group_lost=0
        for i,j in mask:
            tl,tr=ts[i],ts[i+1];ul,ur=us[j],us[j+1]
            rmax=C-F(5,2)*tl-ul;rmin=C-F(5,2)*tr-ur
            if label=='E0':assert rmax/a<=2
            else:assert 2<rmax/a<=F(9,4)
            saved=author.pair_cell(tl,tr,ul,ur,False,False,sf,df,cache)
            assert saved>=0
            if rmin/a<2:
                assert saved==0
                if label=='E1':crossing+=1
            if label=='E0':assert saved==0
            Z=min(ul,F(3,13)*(1-tr-ur));assert a<Z<=ul
            assert 1-tr-ur-F(13,3)*Z>=0
            bracket=Interval(1)/Interval(a)-(log_interval(Z/a)+1)/Interval(Z)
            assert bracket.lo>0
            x=Interval(K)/Interval(80)*lt[i]*lu[j]*bracket
            extra+=x
            inner=lu[j]/Interval(a)+(la[j+1]+2)/Interval(ur)-(la[j]+2)/Interval(ul)
            assert inner.lo>0
            base+=Interval(K)/Interval(8)*lt[i]*inner
            group_lost+=saved
            rows.append([label,i,j,str(Z),str(saved),str(x.lo),str(x.hi)])
        results[label]={'cells':len(mask),'mask_sha256':mask_digest(mask),
            'base_gross_interval':base.record(),'extra_gross_interval':extra.record(),
            'saved_floor_Q80_sum':str(group_lost)}
        all_extra+=extra;lost+=group_lost
    assert crossing==227 and lost==11290806845175652020270
    oldL6=F(assembly['retained_L6_lower_before_new_deletion'])
    floor_loss=F(lost,4*author.Q)
    loss_upper=F(2334885,10**9);old_loss=F(old['adoptable_loss_upper'])
    extra_gross=F(558583,10**9)
    loss_recovery=old_loss-loss_upper
    gain=extra_gross+loss_recovery
    assert floor_loss<loss_upper<old_loss
    assert extra_gross<F(all_extra.lo,P)
    assert gain==F(595792,10**9)>0
    assert oldL6==F(int(read('sharp-248-pair-refinement-512.json')['chosen_lower_numerators']['L6']),author.Q)
    newL6=oldL6-4*loss_upper
    assert newL6>F(assembly['retained_L6_lower'])
    old_gross=F(2095161,10**9)+F(old['adoptable_gross_gain'])
    oldnet=old_gross-old_loss
    newnet=old_gross+extra_gross-loss_upper
    assert newnet-oldnet==gain
    slope_old=F(assembly['feedback_slope']);slope_delta=extra_gross/K
    slope_new=slope_old-slope_delta;assert slope_new>0
    # A direct identity makes duplicate adoption or net/K feedback visible.
    new_raw=F(assembly['raw_c_lower'])+gain
    budget_each=F(assembly['finite_approximation_loss_budget_in_c_units'])+gain/8
    assert F(assembly['finite_approximation_loss_budget_in_c_units'])==F(assembly['compiler_loss_budget_in_c_units'])
    new_after=new_raw-2*budget_each
    assert new_after-F(assembly['raw_after_two_error_budgets'])==3*gain/4
    assert all(sha(n)==digest for n,digest in bindings.items())
    assert state_guard==sha('research-state.json') and main_guard==sha('imbalanced-chen-transfer.tex')
    rowstream=''.join(','.join(map(str,row))+'\n' for row in rows)
    out={'status':'DIRECTED_RATIONAL_SAME_MASK_GROSS_AND_SAVED_LOSS_CERTIFICATE',
      'date':'2026-10-10','adopted':False,'selected_parameters_changed':False,
      'new_cells':0,'old_cells':17572,'historical_total_count_typo':{'printed':17772,'actual':17572},
      'bits':160,'log_working_bits':224,'log_terms':96,
      'K':str(K),'omega_lower':'11/20','omega_argument_threshold':'7/3',
      'inner_fibre':'a<w<z<min(ul,3*(1-tr-ur)/13)',
      'extra_gross_formula':'K/80*log(tr/tl)*log(ur/ul)*(1/a-(1+log(Z/a))/Z)',
      'groups':results,'extra_gross_interval':all_extra.record(),
      'discarded_author_floor_Q80_sum':str(lost),'author_Q':str(author.Q),
      'exact_discarded_E1_floor_in_c_units':str(floor_loss),
      'previous_adopted_loss_upper':str(old_loss),
      'replacement_loss_upper':str(loss_upper),'loss_recovery':str(loss_recovery),
      'extra_gross_lower_adoptable':str(extra_gross),'raw_increment_adoptable':str(gain),
      'previous_total_removed_gross':str(old_gross),'replacement_total_removed_gross':str(old_gross+extra_gross),
      'previous_total_net':str(oldnet),'replacement_total_net':str(newnet),
      'retained_L6_lower_before_crop':str(oldL6),'replacement_retained_L6_lower':str(newL6),
      'feedback_slope_before':str(slope_old),'feedback_slope_deduction_EXTRA_GROSS_over_K':str(slope_delta),
      'feedback_slope_after':str(slope_new),
      'prospective_raw_c_lower':str(new_raw),'prospective_lower_budget_each':str(budget_each),
      'prospective_raw_after_two_error_budgets':str(new_after),
      'extra_allowance_each':str(gain/8),'after_budget_increment':str(3*gain/4),
      'new_c_or_b_selected':False,
      'unchanged_parent_cover':{'E0':{'parents':34405,'stream':cover0['parent_index_stream_sha256']},
        'E1':{'parents':61956,'stream':cover1['parent_index_stream_sha256']},
        'methods':'all switched method0','priority40000_overlap':0},
      'row_columns':['group','t_index','u_index','exact_Z','discarded_floor_Q80','extra_lower_Q160','extra_upper_Q160'],
      'rows':rows,'row_stream_sha256':hashlib.sha256(rowstream.encode()).hexdigest(),
      'literal_author_floor_replay':True,'independent_floor_replay_claimed':False,
      'scope':['New arithmetic proves extra mass and exact saved estimator loss only.',
        'Same current positive-family and restrictable-majorant theorem is retained.',
        'Analytic omega lemma and independent replay have separate proof artifacts.',
        'Old E0/E1 gains are not added again; only the displayed replacement difference is new.',
        'No state, main, old file, formalization or GitHub change.'],
      'source_dependency_pairs_verified':source_maps_checked,
      'bindings_sha256':bindings,'seconds':time.monotonic()-start,
      'display_only':{'extra_gross_lower':float(F(all_extra.lo,P)),
        'exact_floor_loss':float(floor_loss),'new_loss_upper':float(loss_upper),
        'safe_gross_delta':float(extra_gross),'safe_loss_recovery':float(loss_recovery),
        'safe_raw_increment':float(gain),'after_budget_increment':float(3*gain/4),
        'prospective_after_budgets':float(new_after),'new_slope':float(slope_new)}}
    with dest.open('x',encoding='utf-8') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({k:out[k] for k in ['status','old_cells','new_cells',
      'source_dependency_pairs_verified','display_only','seconds']},indent=2),flush=True)

if __name__=='__main__':main()
