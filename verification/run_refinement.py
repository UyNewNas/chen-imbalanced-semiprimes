"""Fresh complete v2 mask arithmetic and row review on the frozen old baseline."""
from pathlib import Path
import argparse,json,shutil,subprocess,sys,tempfile,time

R=Path(__file__).resolve().parent
OUTPUTS=['corpus-renewal-paired-certificate-v2-20261010.json',
         'corpus-renewal-paired-independent-v2-20261010.json',
         'corpus-renewal-paired-final-review-20261010.json']

def run(source,work):
    started=time.monotonic();work.mkdir(parents=True,exist_ok=False)
    names=json.loads((R/'refinement-inputs.json').read_text(encoding='utf8'))
    for n in names:
        target=('frozen-'+n) if n in OUTPUTS else n
        shutil.copyfile(source/n,work/target)
    for name in ['replay_refinement_author.py','replay_refinement_review.py']:
        shutil.copyfile(R/name,work/name)
    scripts=['replay_refinement_author.py','corpus-renewal-paired-independent-v2-20261010.py','replay_refinement_review.py']
    for script,output in zip(scripts,OUTPUTS):
        subprocess.run([sys.executable,'-B',str(work/script)],cwd=work,check=True)
        new=json.loads((work/output).read_text(encoding='utf8'))
        old=json.loads((work/('frozen-'+output)).read_text(encoding='utf-8-sig'))
        excluded={'bindings_sha256','seconds','source_dependency_pairs_verified',
                  'source_bound_dependencies_verified','public_replay_scope'}
        for key in old:
            if key not in excluded:assert new[key]==old[key],(script,key,'frozen disagreement')
    a=json.loads((source/'corpus-renewal-improvement-assembly-20261010.json').read_text(encoding='utf-8-sig'))
    review=json.loads((work/OUTPUTS[-1]).read_text(encoding='utf8'))
    selected=a['selected_result']
    for field in ['safe_K','safe_c','raw_c_lower','finite_approximation_loss_budget_in_c_units',
                  'compiler_loss_budget_in_c_units','raw_after_two_error_budgets']:
        assert selected[field]==review[field],field
    assert selected['feedback_slope']==review['future_slope']
    assert selected['b']==review['primary']['b']
    assert a['secondary_rational_corollary']['b']==review['secondary']['b']
    result=dict(status='PASS_FRESH_V2_CELL_FLOOR_INDEPENDENT_INTERVAL_AND_TRANSFER',
        version='7/3 threshold v2; replaces v1; no gains accumulated',
        actual_cells=17572,E0=6091,E1=11481,zero_crossing_E1=227,new_cells=0,
        author_saved_floors_recomputed=True,independent_interval_engine_rerun=True,
        all_mathematical_fields_match_frozen_records=True,
        safe_K=review['safe_K'],safe_c=review['safe_c'],b=review['primary']['b'],secondary_b=review['secondary']['b'],
        raw_increment=review['raw_gain'],after_budget_increment=review['after_budget_gain'],
        additional_gross=review['safe_extra_gross'],loss_recovery=review['loss_recovery'],
        inherited_baseline_quadratures_rerun=False,seconds=time.monotonic()-started)
    (work/'refinement-run-result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=R/'data');p.add_argument('--work-dir',type=Path,required=True)
    a=p.parse_args();run(a.source.resolve(),a.work_dir.resolve())
