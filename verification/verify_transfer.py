"""Independent exact scalar checks; no floating-point decision is made.

For 0 < b < 1 put t=(1-b)/(1+b). Then
log(1/b)=2 sum_{j>=0} t^(2j+1)/(2j+1).
The omitted positive tail after m terms is at most
2 t^(2m+1)/((2m+1)(1-t^2)).
These checks take K,c and reserve endpoints as inputs; they are not a proof
of the analytic upper and lower sieves that supply those inputs.
"""
from fractions import Fraction as F
from pathlib import Path
import json


def bounds(b, m=96):
    assert 0 < b < 1
    t = (1-b)/(1+b)
    lo = 2*sum((t**(2*j+1)/(2*j+1) for j in range(m)), F(0))
    hi = lo + 2*t**(2*m+1)/((2*m+1)*(1-t*t))
    return lo, hi


def check(source):
    a = json.loads((source/'global-theorem-completion-assembly-20261010.json').read_text(encoding='utf-8-sig'))
    s = a['selected_result']
    K, c = F(s['safe_K']), F(s['safe_c'])
    assert (K,c) == (F('6.7617368'), F('1.96131568'))
    upper_slack = K-F(s['new_upper_raw_endpoint'])-F(s['upper_density_approximation_budget'])-F(s['upper_compiler_budget'])
    lower_slack = F(s['raw_c_lower'])-F(s['finite_approximation_loss_budget_in_c_units'])-F(s['compiler_loss_budget_in_c_units'])-c
    assert upper_slack == F(s['remaining_upper_slack']) > 0
    assert lower_slack == F(s['remaining_lower_slack']) > 0
    rows=[]
    for r, expected_b, expected_count in [(s,F('0.74821909'),F(1,100000)),
        (a['secondary_rational_corollary'],F('0.74821809'),F(1,1000000))]:
        b,count = F(r['b']),F(r['count'])
        assert (b,count)==(expected_b,expected_count)
        lo,hi=bounds(b)
        slack=c-K*hi-count
        assert slack > 0
        prev_lo,prev_hi=bounds(b-F(1,10**8))
        assert c-K*prev_lo-count < 0
        rows.append(dict(b=str(b),count=str(count),positive_slack_exact=str(slack),
                         slack_display_only=float(slack),preceding_eight_decimal_fails=True))
    return dict(status='PASS_EXACT_TRANSFER_AND_RESERVES',K=str(K),c=str(c),
                upper_slack_exact=str(upper_slack),lower_slack_exact=str(lower_slack),
                method='96 terms of positive atanh series with explicit rational tail',
                branches=rows,scope='Scalar implications only; analytic premises are not certified by this script.')


if __name__=='__main__':
    result=check(Path(__file__).resolve().parent/'data')
    print(json.dumps(result,indent=2))
