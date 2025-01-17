# ----------------------------------------------------------------------- #
# How to Lose Some Weight - A Practical Template Syndrome Decoding Attack #
#           code for recreating Tab. 1: the cost estimates                #
# ----------------------------------------------------------------------- #

load('cost_estimates.spyx')

McEliece0 = {'label': 0, 'n': 2197, 'k': 1758, 'w': 37} #our running example
McEliece1 = {'label': 1, 'n': 3488, 'k': 2720, 'w': 64}
McEliece2 = {'label': 2, 'n': 4608, 'k': 3360, 'w': 96}
McEliece3 = {'label': 3, 'n': 6960, 'k': 5413, 'w': 119}
McEliece4 = {'label': 4, 'n': 8192, 'k': 6528, 'w': 128}

paramsets = [McEliece0, McEliece1, McEliece2, McEliece3, McEliece4]
#paramsets = [McEliece0]

Titer  = 0
bitOps = 1
print('~~~~~')
for param in paramsets:
    # bit operations for std Prange #
    n = param['n']
    k = param['k']
    w = param['w']

    ell = 0
    p   = 0

    C_Prange = benchmark_ISD(n, k, w, ell,p, Titer , bitOps)
    print(f'n={n}, C_Prange={C_Prange}')


print('~~~~~')

#optimal parameters, determined via
"""
for param in paramsets:
    # optimizing bit operations for std Dumer #
    n = param['n']
    k = param['k']
    w = param['w']
    
    C_Dumer, p_best, ell_best = optimize_ISD(n, k, w, Titer, bitOps)

    print(f'n = {n}: C_Dumer = {C_Dumer} using ell={ell_best}, p={p_best}')
"""

Dumer0 =  {'label': 0, 'p':  3, 'ell': 29}
Dumer1 =  {'label': 1, 'p':  5, 'ell': 47}
Dumer2 =  {'label': 2, 'p':  5, 'ell': 49}
Dumer3 =  {'label': 3, 'p':  8, 'ell': 78}
Dumer4 =  {'label': 4, 'p': 10, 'ell': 97}
Dumer_params = [Dumer0, Dumer1, Dumer2, Dumer3, Dumer4]

for param in paramsets:
    # bit operations for std Dumer #
    n = param['n']
    k = param['k']
    w = param['w']

    label = param['label']
    Dumer = Dumer_params[label]
    ell = Dumer['ell']
    p = Dumer['p']

    C_Dumer = benchmark_ISD(n, k, w, ell,p, Titer , bitOps)

    print(f'n = {n}: C_Dumer = {C_Dumer} (ell={ell}, p={p})')



bmax = 32
useChecks = 1
advPerm = 1
proportional = 1

niter = 1000 # number of tried templates
print('~~~~~')
for param in paramsets:
    # bit operations for Template Prange
    n = param['n']
    k = param['k']
    w = param['w']

    ell = 0
    p   = 0

    C_TPrange = benchmark_template_ISD(n, k, w, niter, ell, p, bmax, useChecks, advPerm, proportional, Titer, bitOps )
    print(f'n = {n}, C_TPrange = {C_TPrange}')


print('~~~~~')

niter = 1000 # number of tried templates
#optimal parameters, determined via
"""
for param in paramsets:
    # optimizing bit operations for std Dumer #
    n = param['n']
    k = param['k']
    w = param['w']
    
    C_Dumer, p_best, ell_best = optimize_template_ISD(n, k, w, niter,  bmax, useChecks, advPerm, proportional, Titer, bitOps )

    print(f'n = {n}: C_Dumer = {C_Dumer} using ell={ell_best}, p={p_best}')
"""

TDumer0 =  {'label': 0, 'p':  2, 'ell': 16}
TDumer1 =  {'label': 1, 'p':  2, 'ell': 17}
TDumer2 =  {'label': 2, 'p':  3, 'ell': 26}
TDumer3 =  {'label': 3, 'p':  3, 'ell': 28}
TDumer4 =  {'label': 4, 'p':  3, 'ell': 30}
TDumer_params = [TDumer0, TDumer1, TDumer2, TDumer3, TDumer4]

for param in paramsets:
    # bit operations for Template Prange
    n = param['n']
    k = param['k']
    w = param['w']

    label = param['label']
    TDumer = TDumer_params[label]
    ell = TDumer['ell']
    p   = TDumer['p']

    C_TDumer = benchmark_template_ISD(n, k, w, niter, ell, p, bmax, useChecks, advPerm, proportional, Titer, bitOps )
    print(f'n = {n}, C_TDumer = {C_TDumer} (ell={ell}, p={p})')

print('~~~~~')