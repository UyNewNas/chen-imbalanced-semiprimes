"""Audit the full public evidence archive and replay the selected numerical ledger.

This requires the separate full archive. It does not re-run its quadratures.
Only explicitly manifested nonmathematical privacy normalizations are accepted
by the adapted numerical replay. All other SHA-256 comparisons are unchanged.
"""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, tempfile, time
from run_refinement import run as run_refinement


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--work-dir',type=Path,help='New output directory; must not exist.')
    args=parser.parse_args()
    root=args.source.resolve()
    manifest=json.loads((root/'public-full-manifest.json').read_text(encoding='utf-8'))
    t=time.monotonic()
    for row in manifest['files']:
        rel=Path(row['name'])
        assert not rel.is_absolute() and '..' not in rel.parts
        p=root/rel
        assert p.stat().st_size==row['size'],('size',row['name'])
        with p.open('rb') as f:
            actual=hashlib.file_digest(f,'sha256').hexdigest()
        assert actual==row['sha256'],('sha256',row['name'])
    temp=None
    if args.work_dir:
        work=args.work_dir.resolve(); work.mkdir(parents=True,exist_ok=False)
    else:
        temp=tempfile.TemporaryDirectory(prefix='chen-full-verification-')
        work=Path(temp.name)
    # Preserve all archive sources and rerun solely in a fresh copied directory.
    for p in (root/'goldbach-research').rglob('*'):
        if not p.is_file():continue
        target=work/p.relative_to(root/'goldbach-research')
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,target)
    shutil.copyfile(root/'public_normalizations.json',work/'public_normalizations.json')
    shutil.copyfile(root/'replay_public_full.py',work/'replay_public_full.py')
    subprocess.run([sys.executable,'-B',str(work/'replay_public_full.py')],cwd=work,check=True)
    out=json.loads((work/'public-numerical-replay.json').read_text(encoding='utf-8'))
    frozen=json.loads((work/'global-selected-numeric-completion-20261010.json').read_text(encoding='utf-8-sig'))
    # Compare the mathematical output, including all proof-map/count records.
    excluded={'status','bindings_sha256','elapsed_seconds_display_only',
              'public_normalization_exceptions','public_normalization_scope'}
    for key in frozen:
        if key not in excluded:
            assert out[key]==frozen[key],('frozen mathematical result disagreement',key)
    normalized=json.loads((root/'public_normalizations.json').read_text(encoding='utf-8'))['files']
    allowed={r['file'] for r in normalized}
    hits={r['file'] for r in out['public_normalization_exceptions']}
    assert hits.issubset(allowed)
    public_by_name={r['file']:r for r in normalized}
    delta_audit=json.loads((root/'refinement-bindings.json').read_text(encoding='utf8'))
    delta_exceptions=[]
    for row in delta_audit['comparisons']:
        with (root/'goldbach-research'/row['file']).open('rb') as src:
            actual=hashlib.file_digest(src,'sha256').hexdigest()
        if actual!=row['original_sha256']:
            change=public_by_name[row['file']]
            assert (row['original_sha256'],actual)==(change['original_sha256'],change['public_sha256'])
            delta_exceptions.append(dict(record=row['record'],file=row['file']))
    refinement=run_refinement(root/'goldbach-research',work/'refinement')
    result=dict(status='PASS_PUBLIC_FULL_BINDINGS_LEDGER_AND_430_REPLAY',
        manifest_files_checked=len(manifest['files']),
        historical_binding_comparisons=out['direct_dependency_hash_checks'],
        public_normalized_file_count=len(normalized),
        inherited_normalization_comparison_count=len(out['public_normalization_exceptions']),
        inherited_normalization_file_count=len(hits),
        normalization_exceptions=out['public_normalization_exceptions'],
        mathematical_fields_equal_frozen_record=True,
        refinement_binding_comparisons=len(delta_audit['comparisons']),
        refinement_normalization_exceptions=delta_exceptions,
        selected_K=refinement['safe_K'],selected_c=refinement['safe_c'],
        selected_refinement=refinement,
        inherited_baseline_quadratures_rerun=False,seconds=time.monotonic()-t)
    (work/'full-run-result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if args.work_dir: print('Detailed results retained in the requested work directory.')
    if temp:temp.cleanup()


if __name__=='__main__':main()
