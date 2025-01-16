McEliece0 = {'label': 0, 'n': 2197, 'k': 1758, 'w': 37} #our running example
McEliece1 = {'label': 1, 'n': 3488, 'k': 2720, 'w': 64, 'm': 12, 'f':x^12+x^3+1}
McEliece2 = {'label': 2,'n': 4608, 'k': 3360, 'w': 96, 'm': 13, 'f':x^13+x^4+x^3+x+1}
McEliece3 = {'label': 3,'n': 6960, 'k': 5413, 'w': 119,'m': 13, 'f':x^13+x^4+x^3+x+1}
McEliece4 = {'label': 4,'n': 8192, 'k': 6528, 'w': 128,'m': 13, 'f':x^13+x^4+x^3+x+1}


all_mcEliece = [McEliece0, McEliece1, McEliece2, McEliece3, McEliece4]

ns   = [    640,      808,     982,    1101,    1223,    1473,    1665,     1995,     2197,     3488]
ks   = [    512,      647,     786,     881,     979,    1179,    1332,     1596,     1758,     2720]
ws   = [     13,       17,      20,      21,      23,      27,      31,       37,       37,       64]
Ts   = [1704816,  1901275, 2568833, 3405425, 3977964, 5681266, 8143245, 11150513, 11492460, 44893685]
Ts = [log(t*10^(-6)/10000, 2).n() for t in Ts] # log of seconds per iteration

tlength = 32 #template-size (can be 8, 16, 32, 64)
extendH = True #whether we extend H with additional rows
bitCounts = True #adds log(n) to all runtimes

def Prange(paramset):
    """
    Original Prange
    """
    n = paramset['n']
    w = paramset['w']
    k = paramset['k']
    r = optimize_m4ri(n, k)
    RT = log( binomial(n,w),2) - log(binomial(n-k,w) ,2)+log(gauss_complexity(n,k,r),2)
    if  bitCounts: RT+=log(n,2)
    return RT

def optimize_m4ri(n, k):
    """
    helper function for determining the complexity of Gaussian Elemination.
    Borrowed from https://github.com/Crypto-TII/CryptographicEstimators/blob/main/cryptographic_estimators/SDEstimator/sd_helper.py
    """
    mem = oo
    (r, v) = (0, oo)
    for i in range(n - k):
        tmp = log(gauss_complexity(n, k, i), 2)
        if v > tmp and r < mem:
            r = i
            v = tmp
    return r

def gauss_complexity(n, k, r):
    """
    Complexity of Gaussain Elimination
    """
    if r != 0:
        return (r ** 2 + 2 ** r + (n - k - r)) * int(((n + r - 1) / r))
    return (n - k) ** 2

def PuncturedPrange(paramset, wset):
    """
    Complexity of Punctured Prange from
    Vincent Grosso, Pierre-Louis Cayrel, Brice Colombier, and Vlad-Florin Drăgoi.
    "Punctured syndrome decoding problem: Efficient side-channel attacks against classic mceliece"
    """

    n = paramset['n']
    k = paramset['k']
    w = paramset['w']
    m = ceil(paramset['n']/tlength)
    assert(len(wset)==m)
    m0 = 0
    for i in range(m):
        if wset[i]==0: m0+=1
    r = optimize_m4ri(n-m0*tlength, k)
    assert(n>m0*tlength)
    RT = log( binomial(n-m0*tlength,w),2) - log(binomial(n-k,w) ,2)+log(gauss_complexity(n-m0*tlength,k,r),2)
    if bitCounts: RT += log(n,2)
    return RT


def PrangeTemplate(paramset, wset):
    """
    Complexity of Tempate Prange (Algorithm 3)
    """
    n = paramset['n']
    k = paramset['k']
    w = paramset['w']
    codim = n - k
    m = ceil(n/tlength)
    assert(len(wset)==m)
    additional_codim = 0
    if extendH:
        #compute the number of non-empy bins
        for i in range(m):
            if not wset[i]==0: additional_codim+=1
    #print('additional_codim:', additional_codim)

    m0 = m - additional_codim

    codim+=additional_codim
    kset = [0]*m
    rset = [0]*m
    for i in range(m):
        kset[i] = floor(codim*wset[i]/w)
        rset[i] = [(codim*wset[i]/w - kset[i]), i]
    rset_sorted = sorted(rset, key=lambda l:l[0], reverse=True)
    #print(rset_sorted)
    #print(kset)

    assert(sum(kset)<=codim)

    while sum(kset)<codim:
        for el in rset_sorted:
            if sum(kset)==codim: break
            kset[el[1]]+=1
    #print(kset)
    res = 0
    for i in range(m):
        res+=log( binomial(tlength, wset[i]),2) - log( binomial( kset[i] , wset[i]),2)
        #print(wset[i], res.n())


    r = optimize_m4ri(n-m0*tlength, k)
    RT = res+log(gauss_complexity(n-m0*tlength,k,r),2)
    if bitCounts: RT += log(n,2)
    return RT

def get_wset_random(paramset):
    m = ceil(paramset['n']/tlength)
    codim = paramset['n'] - paramset['k']
    wset = [0]*m
    wt = paramset['w']
    while wt > 0:
        pos = ZZ.random_element(0, m)
        assert(wset[pos]<32)
        wset[pos]+=1
        wt -= 1
    assert(sum(wset)==paramset['w'])

    additional_codim = 0
    if extendH:
        #compute the number of non-empy bins
        for i in range(m):
            if not wset[i]==0: additional_codim+=1
    codim+=additional_codim

    #print('additional_codim:', additional_codim)

    wset_scaled = [round(codim*wset[i]/paramset['w']) for i in range(m)]
    return wset

def poisson(mu, r):
    return exp(-mu)*(mu**r) / factorial(r)

def get_wset_ballsAndBins(paramset):
    #due to rounding and the fact the we compute the *expectations*, we usually have sum(wset)!=nballs
    #can correct wset by increasing assigning the remaining balls

    m = ceil(paramset['n']/tlength) #total number of bins
    nballs = paramset['w']
    wset = [0]*m
    full_bins = 0
    i = 0
    additional_codim = 0
    while full_bins < m:
        nbins_with_i_balls = round(m*poisson(nballs/m, i).n())
        for j in range(full_bins, full_bins+nbins_with_i_balls, 1):
            wset[j] = i
        full_bins += nbins_with_i_balls
        #print(i, nbins_with_i_balls, full_bins, m*poisson(nballs/m, i).n())
        if i>0: additional_codim+=nbins_with_i_balls
        i += 1
        if nbins_with_i_balls==0: break
    #print(sum(wset), nballs)

    #print('additional_codim expected:', additional_codim)
    return wset

# Getting Table 1
#"""
for paramset in all_mcEliece:
    avTemplatePrange = 0
    avPuncturedPrange = 0
    Ntrials = 10000
    for i in range(Ntrials):
        wset = get_wset_random(paramset)
        avTemplatePrange+= PrangeTemplate(paramset, wset).n()
        avPuncturedPrange+= PuncturedPrange(paramset, wset).n()
    avTemplatePrange/=Ntrials
    avPuncturedPrange/=Ntrials
    #wset_ = get_wset_ballsAndBins(paramset)
    #T = PrangeTemplate(paramset, wset_).n()

    TPrange = Prange(paramset).n()
    print('Paramset: McEliece', paramset['label'], ' Prange:', TPrange, ' Punctured av.: ', avPuncturedPrange, ' PrangeTemplate av.:', avTemplatePrange)
#"""
"""
for i in range(len(ns)):
    paramset = {'n': ns[i], 'k': ks[i], 'w': ws[i]}
    avTemplatePrange = 0
    avPuncturedPrange = 0
    Ntrials = 5000
    for j in range(Ntrials):
        wset = get_wset_random(paramset)
        avTemplatePrange+= PrangeTemplate(paramset, wset).n()#+Ts[i]
        avPuncturedPrange+=PuncturedPrange(paramset, wset).n()
    avTemplatePrange/=Ntrials
    avPuncturedPrange/=Ntrials
    print('n:', ns[i], 'Prange: ', Prange(paramset).n(), ' PuncturedPrange av:', avPuncturedPrange, ' PrangeTemplate av.:', avTemplatePrange)
"""
