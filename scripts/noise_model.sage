# ----------------------------------------------------------------------- #
# How to Lose Some Weight - A Practical Template Syndrome Decoding Attack #
#           code for recreating Fig. 3: the noise model                   #
# ----------------------------------------------------------------------- #

# running example
n = 2197
w = 37
k = 1758

# CM I
#n = 3488
#w = 64

b = 32

# P(weight w in block)
Pw = [N(binomial(b,ell)*binomial(n-b,w-ell)/binomial(n,w)) for ell in range(b+1)]
print(f'Pw={Pw[:5]}')

# decision boundaries
#sigma = 0.22025
sigma = 0.253
wbar = [N(ell+0.5 + sigma^2*log(Pw[ell]/Pw[ell+1]) ) for ell in range(b)]
print(f'wbar={wbar[:5]}')

import numpy as np
from scipy import special
def qfunc(x):
    return 0.5 - 0.5*special.erf(x/np.sqrt(2))

Perr = Pw[0]*qfunc( (wbar[0]-0)/sigma) + sum(Pw[ell]*(qfunc( (wbar[ell]-ell)/sigma )+qfunc( (ell-wbar[ell-1])/sigma)) for ell in range(1,b)) + Pw[b]*qfunc( (b-wbar[b-1])/sigma)
print(f'Perr={Perr}')
print(f'expected number of errors {Perr*ceil(n/b)}')