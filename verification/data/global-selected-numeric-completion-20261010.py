"""Selected finite numerical proof map, direct ledgers and seed-to-430 replay.

No project Python module is imported. No old artifact or state is modified.
The existing complete quadratures supply mathematical inputs via the explicitly
listed enclosure/coverage proofs; this is not a disguised rerun of quadratures.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib, json, time

B = Path(__file__).resolve().parent
STEM = 'global-selected-numeric-completion-20261010'
Q = 2**80
USED = set()

def read(name):
    USED.add(name)
    return json.loads((B/name).read_text(encoding='utf-8-sig'))

def sha(name):
    with (B/name).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def val(name, key):
    z = read(name)
    for k in key.split('.'):
        z = z[k]
    return F(z)

def log_bounds(x):
    # Independent exact positive Mercator series, no rounded log backend.
    assert 0 < x < 1
    z = 1-x
    lo = sum((z**j/j for j in range(1, 193)), F(0))
    return lo, lo + z**193/(193*x)

def digest(H, h):
    return hashlib.sha256(('H:'+','.join(map(str,H))+'\nh:'+','.join(map(str,h))+'\n').encode()).hexdigest()

def uniform_bank_replay():
    first = read('certify_finite_step_table.json')
    original = read('certify_finite_w2_action_bank.json')
    expanded = read('expanded-w2-midpoint-bank33-v1.json')
    final = read('fixed-action-stationarity-pilot-20261009.json')['arms'][1]
    assert F(final['K']) == F('6.761767')
    constants = final['action_constants']
    actions = []
    for kind, bank in [('W1',first['candidate_bank']),
                       ('W2_original',original['actions']),
                       ('W2_expanded',expanded['actions'])]:
        for index,a in enumerate(bank):
            const = constants[len(actions)]
            assert (const['kind'],const['index']) == (kind,index)
            if kind == 'W1':
                points = [(a['point_index'],1)]
                hw = [(j,int(w)) for j,w in a['feedback_weights']]
                Hw = []
                divisor = 1
                s = F(a['s'])
            else:
                points = a['Gamma1_points']
                hw = [(j,int(lo)) for row in a['negative_h_weights'] for j,lo,hi in row['weights']]
                Hw = [(j,int(lo)) for row in a['ordinary_pair_H_weights'] for j,lo,hi in row['weights']]
                divisor = 5
                s = F(a['parameters'][0])
            assert all(0 <= j < 200 and w >= 0 for j,w in hw+Hw)
            assert a['targets'] == [i for i in range(200) if 1+F(i+1,50) <= s]
            actions.append((int(const['base']),points,hw,Hw,divisor,a['targets']))
    assert len(actions) == 106
    weights = list(map(int,first['log_weight_lower_numerators']))
    assert len(weights) == 200 and min(weights) > 0
    seed = F(2327,100000)
    H = [int(seed*Q) if 1+F(i+1,50) <= F(11,5) else 0 for i in range(200)]
    h = [0]*200
    assert sum(bool(x) for x in H) == 60
    start_hash = digest(H,h)
    comparisons = []
    history = hashlib.sha256()
    for iteration in range(1,431):
        # All right sides read old H,h. Per-product floors are downward.
        sh = [0]*201
        sl = [0]*201
        for j in range(199,-1,-1):
            sh[j] = sh[j+1]+H[j]*weights[j]//Q
            sl[j] = sl[j+1]+h[j]*weights[j]//Q
        high = H.copy()
        low = h.copy()
        for i in range(50,200):
            high[i] = max(high[i], sl[max(50,i-49)])
            low[i] = max(low[i], sh[i-49])
        for base,pts,hw,Hw,divisor,targets in actions:
            z = sum(H[j]*v for j,v in pts)
            z += sum(h[j]*v//Q for j,v in hw)
            z += sum(H[j]*v//Q for j,v in Hw)
            bound = base+z//divisor
            for i in targets:
                high[i] = max(high[i],bound)
        assert all(x >= y >= 0 for x,y in zip(high+low,H+h))
        H,h = high,low
        history.update(f'{iteration}:{digest(H,h)}\n'.encode())
        if iteration in (248,312,430):
            source = ('expanded-w2-midpoint-table248-v1.json' if iteration == 248 else
                      'table-k-feedback-pilot-20261009.json' if iteration == 312 else
                      'fixed-action-stationarity-pilot-20261009.json')
            target = read(source) if iteration != 430 else final
            saved = list(map(int,target['H_cell_lower_numerators']+target['h_cell_lower_numerators']))
            delta = [a-b for a,b in zip(H+h,saved)]
            assert len(delta) == 400 and min(delta) >= 0
            if iteration == 430:
                assert max(delta) == 0 and digest(H,h) == final['fixed_array_sha256']
            comparisons.append(dict(round=iteration,source=source,
                entries=400,shortfalls=0,strict_improvements=sum(d>0 for d in delta),
                all_entries_equal=not any(delta),minimum_integer_difference=str(min(delta)),
                array_sha256=digest(H,h)))
    return dict(K_anchor=str(F(final['K'])),actions=106,rounds=430,
        initial_H_seed=str(seed),initial_nonzero_H_cells=60,initial_h_zero=True,
        input_intermediate_arrays=False,old_transition_module_imported=False,
        original_action_constants_are_proved_numeric_inputs=True,
        comparisons=comparisons,start_sha256=start_hash,
        all430_complete_array_stream_sha256=history.hexdigest(),
        conclusion='The original seed and final finite action bank alone prove the exact selected430 arrays. Earlier selected248/312 arrays are dominated at their respective finite depths. No new arrays or gain are selected.')

def main():
    start = time.monotonic()
    assert not (B/(STEM+'.json')).exists()
    current = read('positive-pair-deletion-improvement-assembly-20261010.json')['selected_result']
    K,c = F(current['safe_K']),F(current['safe_c'])
    assert (K,c) == (F('6.7617368'),F('1.96131568'))
    ext = read('certify_extended_prefix.json')
    upper_components = {
        'ordinary_root':4/F(ext['theta']),
        'first_integral':F(ext['first_integral_upper_exact']),
        'compact_third':F(ext['third_integral_upper_exact']),
        'complete_small_third_tail':val('numerical-primitive-adversarial-20261009.json','independent_tail_upper_Q112')}
    savings = [
      ('replacement_fourth','stronger-fourth-independent-replay.json','gain_lower'),
      ('replacement_old_fifth','independent_fifth_subdivision_replay.json','stored_gain_verified'),
      ('new_disjoint_fifth','fifth-new-coverage-independent-audit.json','gains.new_disjoint'),
      ('replacement_old_sixth','sixth-prefix-expanded-independent-audit.json','gain_lower'),
      ('above_mu_second_coarse_only','second-prefix-upper-fallback-independent-review-20261009.json','no_grid_coarse_lower'),
      ('above_mu_third','third-prefix-upper-fallback-independent-review-20261009.json','selected_lower'),
      ('eleven_sixth_cells_final_log','sixth-log-moment-jensen-independent-20261009.json','new_gain'),
      ('new_descending33','sixth-tent33-descending-independent-20261009.json','gain_lower'),
      ('ascending34_actual_parent','sixth-tent34-ascending-independent-20261009.json','gain_lower'),
      ('descending34_actual_parent','sixth-tent34-descending-independent-20261009.json','gain_lower'),
      ('lower_clip_window','sixth-lower-clip-window-independent-20261009.json','gain_lower'),
      ('lower_descending34_window','sixth-lower-desc34-window-independent-20261009.json','gain_lower'),
      ('selected_local_corrections','sixth-lower-local-additions-independent-20261009.json','selected_gain'),
      ('four_nonsingleton_sixth','sixth-nonsingleton-independent-20261009.json','gain_lower'),
      ('original_eighth_radial_max','eighth-radial-independent-20261009.json','whole_cell_max'),
      ('fourteen_ascending_eighth','eighth-prefix-family-independent-20261009.json','new_gain_lower'),
      ('fifteen_descending_eighth','eighth-prefix-descending-independent-20261009.json','gain_lower'),
      ('fifteen_unequal_eighth_final_log','eighth-unequal-moment-independent-20261009.json','new_gain'),
      ('thirteen_third_minorant_final_log','eighth-third-moment-refinement-independent-20261009.json','new_gain'),
      ('twentynine_history_final_log','eighth-history-moment-independent-20261009.json','new_gain'),
      ('seed_curvature_only_increment','eighth-joint-affine-curvature-independent-20261009.json','increment_over_adopted_log'),
      ('seed_mixed_only_increment','eighth-joint-transverse-mixed-independent-20261009.json','increment')]
    upper_rows = [dict(id=ident,source=name,field=field,value=str(val(name,field))) for ident,name,field in savings]
    E = sum(upper_components.values(),F(0))-sum((F(r['value']) for r in upper_rows),F(0))
    assert E == F(current['new_upper_raw_endpoint'])
    # Every selected upper saving is accounted once, using final values.
    oldledger = read('lichtman-current-upper-density-assembly-20261009.json')
    joint = read('actual-k-density-joint-audit.json')
    Eref = F(joint['endpoint_K_upper'])-F(oldledger['fallbacks']['second'])-F(oldledger['fallbacks']['third'])
    upper_budget = (F('6.7618')-Eref)/4+(Eref-E)/8
    assert upper_budget == F(current['upper_density_approximation_budget']) == F(current['upper_compiler_budget'])
    upper_slack = K-E-2*upper_budget
    assert upper_slack == F(current['remaining_upper_slack']) > 0
    compact = read('numerical-third-compact-independent-full-20261010.json')
    assert F(compact['upper_exact']) <= F(ext['third_integral_upper_exact'])

    # Direct eleven-term enclosure: no old n32 cost or old four-prefix patch.
    baseline = read('sharp-248-baseline-128.json')['terms']
    def endpoint(name,side):
        z = baseline[name]
        return F(int(z[side+'_numerator']),int(z['denominator']))
    pairs = read('independent-sharp-248-pair512-replay-20261008.json')
    assert pairs['cell_counts'] == {'L5':512*513//2,'L6':2*512*512}
    switched = read('sharp-switch-midpoint-1024.json')['terms']
    def switch(name):
        z = switched[name]
        return F(int(z['upper_numerator']),int(z['denominator']))
    n48 = read('sharp-248-four-min-48.json')
    n48check = read('independent-sharp-248-four-min48-replay-20261008.json')
    assert n48check['all_rows_checked'] == comb(51,4)+48*comb(50,3) == 1190700
    assert not n48check['floor_shortfalls']
    terms = dict(L1=endpoint('L1','lower'),L2=endpoint('L2','lower'),
        L5=F(pairs['selected_lower']['L5']),L6=F(pairs['selected_lower']['L6']),
        U3=endpoint('U3','upper'),U4=endpoint('U4','upper'),
        U7=switch('U7'),U8=switch('U8'),U9=switch('U9'),
        U10_plus_U11=F(n48['U10_U11_cost_upper']))
    weights = dict(L1=3,L2=1,L5=1,L6=1,U3=-1,U4=-1,U7=-2,U8=-1,U9=-1,U10_plus_U11=-1)
    eleven = sum((weights[n]*v for n,v in terms.items()),F(0))/4
    priority = read('four-min-priority40000-independent-20261009.json')
    assert priority['all_parents'] == 40000 and priority['all_children'] == 516779
    assert priority['child_floor_shortfalls'] == 0
    anchor = read('lower-start-convex-corner-transfer-independent-review.json')
    ordinary = F(anchor['independent_ordinary_interval'][0])
    removed = F(anchor['independent_removed_H248_interval'][1])
    lower = read('lower-start-method-coverage-20261009.json')['selected_result']
    lower_integrals = {n:F(lower[n]) for n in ('I3_a','I3_h','I4_a','I4_h','I5_a','I5_h')}
    Iweighted = sum((v*(3 if n.endswith('_a') else 1) for n,v in lower_integrals.items()),F(0))
    restoration = val('lower-fourth-restoration-slanted-independent-20261009.json','author_lower_covered')
    base = eleven+F(priority['coefficient_gain_lower'])+(ordinary-removed+Iweighted+restoration)/4
    base_selected = val('slanted-restoration-assembly-20261009.json','selected_result.raw_c_lower')
    assert base >= base_selected
    slope = read('selected-sieve-budget-feedback-audit-20261010.json')['slope']
    A = F(slope['total'])
    assert A == F(slope['U789'])+F(slope['retained_U1011'])
    t312 = read('table-k-feedback-transfer-assembly-20261009.json')['selected_result']
    gain312 = F(t312['U3_U4_table_gain'])
    gain430 = val('fixed-action-current-k-two-cost-20261009.json','gain_in_c_lower')
    zero_gain = F('0.002095161')
    gross,loss = F('0.003799514'),F('0.002372094')
    crop = read('positive-pair-deletion-independent-rational-20261010.json')
    assert crop['positive_count'] == 11481 and crop['unresolved_count'] == 0
    assert F(crop['gain_interval']['lower']) > gross
    assert F(crop['old_floor_loss_upper_proxy_interval']['upper']) < loss
    direct_raw = base+A*(F('6.7619')-K)+gain312+gain430+zero_gain+gross-loss
    assert direct_raw-F(current['raw_c_lower']) == base-base_selected >= 0
    # Old budgets are retained; no gain is used to hide an enclosure shortfall.
    lower_budget = F(t312['finite_approximation_loss_budget_in_c_units'])+(A*(F('6.7618')-K)+gain430+zero_gain+gross-loss)/8
    assert lower_budget == F(current['finite_approximation_loss_budget_in_c_units']) == F(current['compiler_loss_budget_in_c_units'])
    lower_slack = F(current['raw_c_lower'])-2*lower_budget-c
    assert lower_slack == F(current['remaining_lower_slack']) > 0
    assert F(current['retained_L6_lower']) == terms['L6']-4*loss
    assert F(current['feedback_slope']) == A-(zero_gain+gross)/K
    transfer = []
    for b,count in [(F('0.74821909'),F(1,100000)),(F('0.74821809'),F(1,1000000))]:
        lo,hi = log_bounds(b)
        surplus = c-K*hi-count
        assert surplus > 0
        transfer.append(dict(b=str(b),count=str(count),positive_surplus_exact=str(surplus),
            surplus_display_only=float(surplus)))
    seed = read('replay_second_seed_first_feedback.json')
    assert F(seed['ordinary_baseline']['lower']) > F(20249,10**6) > F(1,50)
    assert F(seed['feedback_coefficient_lower']) > F(13,100)
    J = F(1,50)
    for _ in range(4): J = F(20249,10**6)+F(13,100)*J
    assert J > F(2327,100000)
    recurrence = uniform_bank_replay()

    # These are proof/result pairs, not an assertion that every history file is
    # needed. Names omitted from this map are not additional numeric premises.
    proof_map = read(STEM+'-obligations.json')
    for item in proof_map['obligations']:
        assert item['status'] == 'covered_by_proved_enclosure_and_complete_selected_output'
        for name in item['proof_sources']+item['complete_records']:
            assert (B/name).is_file(),name
            USED.add(name)
    # Confirm actual output coverage, rather than status labels alone.
    counts = {
      'upper_fourth_boxes':read('stronger-fourth-independent-replay.json')['boxes'],
      'upper_fifth_old_parents':read('independent_fifth_subdivision_replay.json')['verified_parents'],
      'upper_fifth_old_children':read('independent_fifth_subdivision_replay.json')['verified_children'],
      'upper_fifth_new_boxes':read('fifth-new-coverage-independent-audit.json')['replay']['all_cells'],
      'upper_old_sixth_boxes':read('sixth-prefix-expanded-independent-audit.json')['cells_verified'],
      'W1_actions':read('numerical-w1-primitive-independent-20261010.json')['counts']['W1_actions'],
      'W2_actions':read('numerical-w2-affine-series-independent-20261010.json')['action_count'],
      'ordinary_function_node_enclosures':read('numerical-lower-shared-sieve-chord-independent-20261011.json')['total_function_nodes_checked'],
      'lower_restoration_boxes':read('lower-fourth-restoration-slanted-independent-20261009.json')['boxes_checked'],
      'lower_restoration_layers':read('lower-fourth-restoration-slanted-independent-20261009.json')['layers_checked']}
    assert counts == dict(upper_fourth_boxes=3699127,upper_fifth_old_parents=292633,
        upper_fifth_old_children=18728512,upper_fifth_new_boxes=1763420,upper_old_sixth_boxes=632104,
        W1_actions=70,W2_actions=36,ordinary_function_node_enclosures=23199750,
        lower_restoration_boxes=98668,lower_restoration_layers=789344)
    assert counts['upper_fifth_old_children'] == 64*counts['upper_fifth_old_parents']
    assert read('third-m1600-arrangement-independent-20261010.json')['complete_both_domains']
    assert read('third-m1600-arrangement-independent-20261010.json')['limit'] == 0
    assert sha('imbalanced-chen-transfer.tex') == '06ecf0c7a8d767115dd0ef595e52844d49486b4b98d5979300faf47f4e09ff78'
    # Bind each supplied direct dependency. Hashes identify evidence, not proofs.
    hashes = {n:sha(n) for n in sorted(USED)}
    direct_binding_checks = 0
    for name in sorted(USED):
        if not name.endswith('.json'):continue
        obj = json.loads((B/name).read_text(encoding='utf-8-sig'))
        for key in ('bindings','bindings_sha256','dependencies_sha256','verified_input_hashes','sources_sha256'):
            refs = obj.get(key,{})
            if not isinstance(refs,dict):continue
            for dep,digest0 in refs.items():
                if not isinstance(digest0,str) or len(digest0)!=64:continue
                path = B/dep
                if not path.is_file():
                    # These selected numeric sources all bind project-relative files.
                    raise AssertionError((name,dep,'missing dependency'))
                actual = hashes.setdefault(dep,sha(dep))
                assert actual == digest0,(name,dep)
                direct_binding_checks += 1
    hashes[Path(__file__).name] = sha(Path(__file__).name)
    hashes[STEM+'.txt'] = sha(STEM+'.txt')
    out = dict(status='SELECTED_FINITE_NUMERICAL_OBLIGATIONS_COVERED_WITH_DIRECT_LEDGER_AND_SEED_TO_430_REPLAY',
      date='2026-10-10',selected_K=str(K),selected_c=str(c),
      scope='Finite numerical inequalities and selected integration/selection coverage only; named analytic function identities, physical density interpretation, legal coefficient construction, PP and Chen weight are separate proof requirements.',
      selected_upper_direct_base={n:str(v) for n,v in upper_components.items()},
      selected_upper_direct_savings=upper_rows,upper_saving_term_count=len(upper_rows),
      upper_raw_exact=str(E),upper_budget_each=str(upper_budget),upper_slack=str(upper_slack),
      direct_lower_terms={n:dict(endpoint=str(v),weight=weights[n]) for n,v in terms.items()},
      direct_lower_initial_eleven=str(eleven),ordinary_start_change_lower=str(ordinary),removed_H248_upper=str(removed),
      lower_selected_integrals={n:str(v) for n,v in lower_integrals.items()},restoration_L2_lower=str(restoration),
      direct_raw_at_K0=str(base),conservative_raw_at_K0=str(base_selected),
      feedback_original_slope=str(A),table_gain248_to312=str(gain312),table_gain312_to430=str(gain430),
      direct_current_raw=str(direct_raw),selected_current_raw=current['raw_c_lower'],
      retained_rounding_surplus=str(direct_raw-F(current['raw_c_lower'])),
      lower_budget_each=str(lower_budget),lower_slack=str(lower_slack),
      transfer=transfer,seed_to_final_bank=recurrence,complete_count_controls=counts,
      obligations=proof_map['obligations'],unmet_selected_finite_numerical_items=[],
      history_only_or_unselected_excluded=proof_map['excluded'],
      fresh_quadrature_rerun=False,new_gain_or_constants_selected=False,
      analytic_requirements_not_certified_here=proof_map['analytic_outside_scope'],
      direct_dependency_hash_checks=direct_binding_checks,
      bindings_sha256=hashes,elapsed_seconds_display_only=time.monotonic()-start)
    (B/(STEM+'.json')).write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=out['status'],obligations=len(out['obligations']),
        upper_terms=len(upper_rows),direct_lower_surplus=str(direct_raw-F(current['raw_c_lower'])),
        seed_to430=recurrence['comparisons'],bindings=len(hashes),seconds=out['elapsed_seconds_display_only']),indent=2))

if __name__ == '__main__':main()
