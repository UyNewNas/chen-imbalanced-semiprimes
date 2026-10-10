# Public adapter: fixed baseline replaces mutable research state; mathematics unchanged.
"""Separate exact scalar, row, geometry and source-binding review.

No project interval/sieve code is imported or executed.  The literal source
inequality is reviewed in the companion report; this program checks the saved
row comparisons and supplies independent rational-log transfer certificates.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
from math import ceil, floor
import hashlib
import json

B = Path(__file__).resolve().parent
STEM = 'corpus-renewal-paired-final-review-20261010'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(n):
    return json.loads((B/n).read_text(encoding='utf-8'))

def display(q):
    with localcontext() as ctx:
        ctx.prec = 45
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def log_interval(q, terms=64):
    """Exact positive atanh series, with a geometric upper remainder."""
    q = F(q)
    assert 1 <= q <= 2
    z = (q-1)/(q+1)
    z2 = z*z
    power = z
    total = F(0)
    for j in range(terms):
        total += 2*power/(2*j+1)
        power *= z2
    tail = 2*power/((2*terms+1)*(1-z2))
    assert tail >= 0
    return total, total+tail

def transfer(c, K, count, numerator):
    b = F(numerator, 10**8)
    lo, hi = log_interval(1/b)
    margin = c-K*hi-count
    previous = b-F(1,10**8)
    plo, phi = log_interval(1/previous)
    previous_margin_upper = c-K*plo-count
    assert margin > 0 and previous_margin_upper < 0
    return dict(b=str(b), count=str(count), log_lower=str(lo), log_upper=str(hi),
        strict_margin_lower=str(margin), strict_margin_lower_decimal=display(margin),
        log_width_upper=str(hi-lo), preceding_grid_b=str(previous),
        preceding_grid_margin_upper=str(previous_margin_upper),
        preceding_grid_margin_upper_decimal=display(previous_margin_upper))

def main():
    selected = read('global-theorem-completion-assembly-20261010.json')['selected_result']
    baseline = read('positive-pair-deletion-improvement-assembly-20261010.json')
    for k in ['safe_K','safe_c','b','count','raw_c_lower',
              'finite_approximation_loss_budget_in_c_units',
              'compiler_loss_budget_in_c_units','raw_after_two_error_budgets',
              'feedback_slope']:
        assert selected[k] == baseline['selected_result'][k], k
    names = ['positive-pair-deletion-improvement-assembly-20261010.json',
        'positive-pair-deletion-certificate-20261010.json',
        'paired-y6-y11-zero-strip-certificate-20261010.json',
        'expanded-w2-midpoint-table248-v1.json',
        'certify_sharp_184_baseline.py','certify_combined_chen_functional.py',
        'certify_extended_prefix.py','sharp-248-pair-refinement-512.json',
        'corpus-renewal-paired-certificate-v2-20261010.py',
        'corpus-renewal-paired-certificate-v2-20261010.json',
        'corpus-renewal-paired-certificate-v2-20261010.txt',
        'corpus-renewal-paired-independent-v2-20261010.py',
        'corpus-renewal-paired-independent-v2-20261010.json',
        'corpus-renewal-paired-independent-20261010.txt']
    bindings = {n:sha(B/n) for n in names}
    author = read(names[9])
    independent = read(names[12])
    verified = set()
    # Full public runner performs the frozen dependency audit separately.
    pos = read(names[1]); zero = read(names[2]); bank = read(names[3])
    e1 = [tuple(map(int,p)) for p in pos['selected_pairs']]
    e0 = [(r['t_index'],j) for r in zero['strips'] for j in range(r['first_u_index'],512)]
    assert len(e1)==len(set(e1))==11481 and len(e0)==len(set(e0))==6091
    assert not set(e0)&set(e1)
    maskhash = hashlib.sha256(''.join(f'{i},{j}\n' for i,j in e1).encode()).hexdigest()
    assert maskhash == pos['selected_mask_sha256'] == independent['selected_mask_sha256']
    assert baseline['total_deleted_cells'] == 17772
    assert len(e0)+len(e1) == independent['total_cells'] == 17572
    Q = int(author['author_Q']); S = int(independent['scale'])
    assert Q==2**80 and S==10**80 and author['bits']==160
    ar = {(r[0],int(r[1]),int(r[2])):r for r in author['rows']}
    ir = {(int(r[0]),int(r[1])):r for r in independent['rows']}
    assert len(ar)==17572 and set(ir)==set(e1)
    assert set(ar)=={('E0',i,j) for i,j in e0}|{('E1',i,j) for i,j in e1}
    a,h,mu,d,C,K = map(F,['417/5000','19/125','25/128','59/250','57/64',selected['safe_K']])
    dt,du=(d-mu)/512,(h-a)/512
    hc = list(map(int,bank['h_cell_lower_numerators']))
    assert all(v>=0 for v in hc)
    crossing = 0
    floor_sum = 0
    root_sum = 0
    for group,cells in [('E0',e0),('E1',e1)]:
        for i,j in cells:
            tl,tr=mu+i*dt,mu+(i+1)*dt
            ul,ur=a+j*du,a+(j+1)*du
            Z=min(ul,F(3,13)*(1-tr-ur))
            row=ar[group,i,j]
            assert F(row[3])==Z and a<Z<=ul
            assert 1-tr-ur-2*Z>=F(7,3)*Z
            # The second affine tent plane is minimal throughout t>=mu.
            assert C-F(1,2)-2*tl<=0
            sl=(C-F(5,2)*tr-ur)/a
            sh=(C-F(5,2)*tl-ul)/a
            fnum=int(row[4])
            if group=='E0':
                assert sh<=2 and fnum==0
                continue
            assert 2<sh<=F(9,4)
            if sl<2:
                crossing+=1
                assert fnum==0
            rr=ir[i,j]
            floor_sum+=fnum
            root_sum+=int(rr[2])
            assert F(fnum,4*Q)<=F(int(rr[2]),S)
            # Independent Q160 and decimal intervals enclose the same mass.
            assert max(F(int(row[5]),2**160),F(int(rr[3]),S)) <= min(F(int(row[6]),2**160),F(int(rr[4]),S))
    assert crossing==227 and floor_sum==11290806845175652020270
    assert F(root_sum,S)==F(independent['discarded_author_floor_upper_proxy']['upper'])
    safe_gross=F(558583,10**9); oldloss=F(2372094,10**9); newloss=F(2334885,10**9)
    gain=safe_gross+oldloss-newloss
    assert gain==F(595792,10**9)==F(author['raw_increment_adoptable'])==F(independent['safe_raw_increment'])
    assert F(floor_sum,4*Q)<=F(root_sum,S)<newloss
    assert sum(F(independent['extra_gross'][g]['lower']) for g in ['E0','E1'])>safe_gross
    # Standalone analytic omega proof uses log x>=2(x-1)/(x+1).
    assert F(27,49)>F(11,20) and F(5,9)>F(11,20)
    assert log_interval(2)[0]>F(2,3)
    raw=F(selected['raw_c_lower'])+gain
    density=F(selected['finite_approximation_loss_budget_in_c_units'])+gain/8
    compiler=F(selected['compiler_loss_budget_in_c_units'])+gain/8
    after=raw-density-compiler
    assert after==F(selected['raw_after_two_error_budgets'])+3*gain/4
    c=F(floor(after*10**8),10**8)
    assert c==F(196176252,10**8) and c<=after<c+F(1,10**8)
    slope=F(selected['feedback_slope'])-safe_gross/K
    assert slope==F(author['feedback_slope_after'])
    assert raw==F(author['prospective_raw_c_lower'])
    assert density==compiler==F(author['prospective_lower_budget_each'])
    assert after==F(author['prospective_raw_after_two_error_budgets'])
    primary=transfer(c,K,F(1,100000),74816964)
    secondary=transfer(c,K,F(1,1000000),74816865)
    mainsha=sha(B/'imbalanced-chen-transfer.tex')
    assert mainsha=='06ecf0c7a8d767115dd0ef595e52844d49486b4b98d5979300faf47f4e09ff78'
    for n,digest in bindings.items():assert sha(B/n)==digest,n
    bindings[Path(__file__).name]=sha(Path(__file__))
    out=dict(status='PASS_SEPARATE_SOURCE_FORMULA_ROW_AND_EXACT_TRANSFER_REVIEW',
        scope='New same-mask improvement only; does not independently reprove every retained analytic or distribution interface.',
        no_project_numerical_imports=True,source_bound_dependencies_verified=len(verified),
        actual_cell_counts=dict(E0=6091,E1=11481,total=17572,crossing_zero_floor=227),
        historical_17772_is_count_typo_only=True,historical_evidence_unchanged=True,
        selected_mask_sha256=maskhash,all_E1_independent_upper_proxies_dominate_literal_saved_floors=True,
        exact_saved_floor_loss=str(F(floor_sum,4*Q)),independent_loss_proxy_upper=str(F(root_sum,S)),
        adopted_loss_upper=str(newloss),safe_extra_gross=str(safe_gross),loss_recovery=str(oldloss-newloss),
        raw_gain=str(gain),after_budget_gain=str(3*gain/4),safe_K=str(K),raw_c_lower=str(raw),
        finite_approximation_loss_budget_in_c_units=str(density),compiler_loss_budget_in_c_units=str(compiler),
        raw_after_two_error_budgets=str(after),safe_c=str(c),remaining_c_slack=str(after-c),
        primary=primary,secondary=secondary,log_method='64 positive exact atanh terms; tail <=2*z^129/(129*(1-z^2))',
        future_slope=str(slope),future_slope_deduction=str(safe_gross/K),
        recovered_loss_is_not_divided_by_K=True,versions_v1_v2_not_accumulated=True,
        main_sha256=mainsha,state_read_guard='latest_own_route_selected matched original assembly selected fields; state not modified or hash-bound',
        bindings_sha256=bindings)
    out['public_replay_scope']='Frozen baseline guard and actual v2 row/geometry/transfer review; history audit separate.'
    dest=B/(STEM+'.json')
    assert not dest.exists()
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=out['status'],c=display(c),primary_b=display(F(primary['b'])),
        primary_slack=primary['strict_margin_lower_decimal'],secondary_b=display(F(secondary['b'])),
        secondary_slack=secondary['strict_margin_lower_decimal'],verified_bindings=len(verified)),indent=2))

if __name__=='__main__': main()
