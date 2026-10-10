"""Validate immutable compact inputs and replay into a fresh private directory."""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, tempfile, time
from verify_latest_transfer import check
from run_refinement import run as run_refinement

ROOT=Path(__file__).resolve().parent


def validate(root, manifest_name='manifest.json'):
    manifest=json.loads((root/manifest_name).read_text(encoding='utf-8'))
    for row in manifest['files']:
        relative=Path(row['name'])
        assert not relative.is_absolute() and '..' not in relative.parts
        p=root/relative
        assert p.stat().st_size==row['size'], ('size',row['name'])
        with p.open('rb') as f:
            actual=hashlib.file_digest(f,'sha256').hexdigest()
        assert actual==row['sha256'],('sha256',row['name'])
    return len(manifest['files'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir','--work',type=Path,help='New directory to retain detailed results; must not exist.')
    args=parser.parse_args()
    started=time.monotonic()
    checked=validate(ROOT)
    temp=None
    if args.work_dir:
        work=args.work_dir.resolve()
        work.mkdir(parents=True,exist_ok=False)
    else:
        temp=tempfile.TemporaryDirectory(prefix='chen-verification-')
        work=Path(temp.name)
    names=json.loads((ROOT/'quick-inputs.json').read_text(encoding='utf-8'))
    for n in names:
        shutil.copyfile(ROOT/'data'/n,work/n)
    shutil.copyfile(ROOT/'data'/'global-selected-numeric-completion-20261010.json',work/'frozen-numerical-result.json')
    shutil.copyfile(ROOT/'replay_selected_ledger.py',work/'replay_selected_ledger.py')
    scalar=check(ROOT/'data')
    (work/'transfer-replay.json').write_text(json.dumps(scalar,indent=2)+'\n',encoding='utf-8')
    subprocess.run([sys.executable,'-B',str(work/'replay_selected_ledger.py')],cwd=work,check=True)
    ledger=json.loads((work/'ledger-replay.json').read_text(encoding='utf-8'))
    refinement=run_refinement(ROOT/'data',work/'refinement')
    result=dict(status='PASS_COMPACT_VERIFICATION',manifest_files_checked=checked,
        ledger_input_files=len(names),upper_savings=ledger['upper_saving_term_count'],
        actions=ledger['seed_to_final_bank']['actions'],rounds=ledger['seed_to_final_bank']['rounds'],
        final_array_sha256=ledger['seed_to_final_bank']['comparisons'][-1]['array_sha256'],
        scalar_status=scalar['status'],
        transfer=[{k:r[k] for k in ('b','count','slack_display_only','preceding_eight_decimal_fails')} for r in scalar['branches']],
        selected_K=refinement['safe_K'],selected_c=refinement['safe_c'],
        selected_refinement=refinement,
        inherited_baseline_quadratures_rerun=False,full_historical_binding_audit_performed=False,
        seconds=time.monotonic()-started)
    (work/'run-result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if args.work_dir: print('Detailed results retained in the requested work directory.')
    if temp: temp.cleanup()


if __name__=='__main__': main()
