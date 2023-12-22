import math
from math import log2, inf
import random
from Instance import * # everything, now an instance can be accessed just by typing inst156
#import Instance # just imports the Instance module, a variable can be accessed with Instance.inst156

# just effectively renaming the comb function
def binom(n: int, k: int):
    return math.comb(n, k)

# -------------------------- NECCESSARY HELPER STUFF ------------------------

# A Block containing columns, necessary for complexity analysis of actual instances.
# Same structure as in the c++ code.
class ColumnsBlock:
    def __init__(self, startpos, length, weight, nrCols = 0):
        self.length = length
        self.startpos = startpos
        self.weight = weight
        self.nrCols = nrCols

    def print(self):
        print("{} {} {} {}".format(self.startpos, self.length, self.nrCols, self.weight))

# Takes in an actual ReducedInstance from Instance.py and parses the error weight string in the same manner
# as the c++ code. Returns a list of blocks
def parseWeightString(inst: ReducedInstance):
    weightString = inst.weightDistro
    ret = []
    for i in range(len(weightString)-1):
        if weightString[i] != "0":
            ret.append(ColumnsBlock(i*32, 32, int(weightString[i])))
    
    m = 32 if inst.n%32 == 0 else inst.n%32

    if weightString[-1:] != "0":
        ret.append(ColumnsBlock((len(weightString)-1)*32, m, int(weightString[-1:])))

    # for b in ret:
    #     b.print()

    # Fancy new computation of nrColsToChose
    tupleList = []
    nkm = inst.n - inst.k + inst.addRows
    wD = inst.w
    sum = 0
    for block in ret:
        resF = (block.weight / wD) * nkm
        res = math.floor(resF)
        block.nrCols = res
        sum += res
        tupleList.append([resF-res, block])
    
    # for b in ret:
    #     b.print()

    tupleList.sort(key=lambda x: x[0], reverse=True)

    while sum < nkm:
        for t in tupleList:
            block = t[1]
            if block.nrCols < block.length:
                block.nrCols = block.nrCols+1
                sum += 1
            if sum == nkm:
                break

    # for b in ret:
    #     b.print()    

    return ret


# ------------------------ PRANGE COMPLEXITY ------------------------

# Computes the number of expected permutations for prange using a standard random permutation.
def prangeComplexityRandomPerm(cols, rows, weight, log = False):
    combs = math.comb(cols, weight) / math.comb(rows, weight)
    ret = math.log2(combs) if log else combs
    return ret

# Same as above, but takes in an Instance object.
def prangeComplexityRandomPerm(inst: Instance, log = True):
    combs = math.comb(inst.ncols, inst.w) / math.comb(inst.nrows, inst.w)
    ret = math.log2(combs) if log else combs
    return ret

# Computes the number of expected permutations for prange for a given actual ReducedInstance.
# Uses a better permutation by taking the weights of the blocks into consideration. 
def prangeComplexitySpecialPerm(inst: ReducedInstance, log = True):
    blocks = parseWeightString(inst)
    exp_perm = 1.0
    for block in blocks:
        t = math.comb(block.nrCols, block.weight)
        b = math.comb(block.length, block.weight)
        exp_perm *= (t/b)
    exp_perm = 1.0 / exp_perm
    exp_perm = log2(exp_perm) if log else math.ceil(exp_perm)
    #print(math.ceil(exp_perm))
    return exp_perm


# ------------------------- DUMER COMPLEXITY ------------------------
# Most of these function helper functions.

def _gaussian_elimination_complexity(n, k, r):
    if r != 0:
        return (r ** 2 + 2 ** r + (n - k - r)) * int(((n + r - 1) / r))

    return (n - k) ** 2

def _gaussian_elimination_complexity2(ncols, nrows, r):
    if r != 0:
        return (r ** 2 + 2 ** r + (nrows - r)) * int(((ncols + r - 1) / r))

    return (nrows) ** 2

def _optimize_m4ri(n, k, ar=0, mem=inf):
    (r, v) = (0, inf)
    for i in range(n - k):
        tmp = math.log2(_gaussian_elimination_complexity(n, k, i))
        if v > tmp and r < mem:
            r = i
            v = tmp
    return r

def _optimize_m4ri2(ncols, nrows, mem=inf):
    (r, v) = (0, inf)
    for i in range(nrows):
        tmp = math.log2(_gaussian_elimination_complexity(ncols, nrows, i))
        if v > tmp and r < mem:
            r = i
            v = tmp
    return r

def _list_merge_complexity(L, l, hmap):
    """
    Complexity estimate of merging two lists exact
    INPUT:
    - ``L`` -- size of lists to be merged
    - ``l`` -- amount of bits used for matching
    - ``hmap`` -- indicates if hashmap is being used (Default 0: no hashmap)
    """
    if L == 1:
        return 1
    if not hmap:
        return max(1, 2 * int(log2(L)) * L + L ** 2 // 2 ** l)
    else:
        return 2 * L + L ** 2 // 2 ** l

def __memory_access_cost(mem, memory_access):
    if memory_access == 0:
        return 0
    elif memory_access == 1:
        return log2(mem)
    elif memory_access == 2:
        return mem / 2
    elif memory_access == 3:
        return mem / 3
    elif callable(memory_access):
        return memory_access(mem)
    return 0

def _mem_matrix(n, k, r):
    """
    Memory usage of parity check matrix in vector space elements
    INPUT:
    - ``n`` -- length of the code
    - ``k`` -- dimension of the code
    - ``r`` -- block size of M4RI procedure
    """
    return n - k + 2 ** r

def _mem_matrix2(nrows, r):
    """
    Memory usage of parity check matrix in vector space elements
    INPUT:
    - ``n`` -- length of the code
    - ``k`` -- dimension of the code
    - ``r`` -- block size of M4RI procedure
    """
    return nrows + 2 ** r

# More or less reference function that was modified. Use function below.
def dumer_complexity(n: int, k: int, w: int, mem=inf, memory_access=0, hmap=1,
                     val_l=0, val_p=0):
    solutions = max(0., log2(binom(n, w)) - (n - k))
    time = inf
    memory = 0
    r = _optimize_m4ri(n, k, mem)

    i_val = [10, 40]
    i_val_inc = [10, 10]
    params = [-1 for _ in range(2)]
    while True:
        stop = True
        for p in range(min(w // 2, i_val[0])):
            # if a p value is given, make sure to set it
            if val_p and p != val_p:
                continue

            for l in range(min(n - k - (w - p), i_val[1])):
                if val_l and l != val_l:
                    continue

                k1 = (k + l) // 2
                L1 = binom(k1, p)
                if log2(L1) > time:
                    continue

                tmp_mem = log2(2 * L1 + _mem_matrix(n, k, r))
                if tmp_mem > mem:
                    continue
                else:
                    Tp = max(log2(binom(n, w)) - log2(binom(n-k-l, w - 2 * p))
                             - log2(binom(k1, p) ** 2) - solutions, 0.)

                #Tp = max(log2(binom(n, w)) - log2(binom(n - k, w - 2 * p)) - log2(binom(k1, p) ** 2) - solutions, 0.)
                Tg = _gaussian_elimination_complexity(n, k, r)
                tmp = Tp + log2(Tg + _list_merge_complexity(L1, l, hmap))

                tmp += __memory_access_cost(tmp_mem, memory_access)

                time = min(time, tmp)
                if tmp == time:
                    memory = tmp_mem
                    params = [p, l, Tp]

        for i in range(len(i_val)):
            if params[i] == i_val[i] - 1:
                stop = False
                i_val[i] += i_val_inc[i]

        if stop:
            break

    par = {"l": params[1], "p": params[0], "c": 0}
    res = {"time": time, "perms": params[2], "memory": memory, "parameters": par}
    #print(res)
    return res

# This function calculates the "optimal" l and p for a given instance when applying
# a standard permutation. Works for both a normal and a reduced instance.
def dumer_complexity2(I: Instance, mem=inf, memory_access=0, hmap=1,
                     val_l=0, val_p=0):
    ncols = I.ncols
    nrows = I.nrows
    w = I.w

    solutions = max(0., log2(binom(ncols, w)) - (nrows))
    time = inf
    memory = 0
    r = _optimize_m4ri2(ncols, nrows, mem)

    i_val = [10, 40]
    i_val_inc = [10, 10]
    params = [-1 for _ in range(2)]
    while True:
        stop = True
        for p in range(min(w // 2, i_val[0])):
            # if a p value is given, make sure to set it
            if val_p and p != val_p:
                continue

            for l in range(min(nrows - (w - p), i_val[1])):
                if val_l and l != val_l:
                    continue

                right_side = I.ncols - (I.nrows -l )
                k1 = right_side // 2
                L1 = binom(k1, p)
                if log2(L1) > time:
                    continue

                tmp_mem = log2(2 * L1 + _mem_matrix2(nrows, r))
                if tmp_mem > mem:
                    continue
                else:
                    Tp = max(log2(binom(ncols, w)) - log2(binom(nrows-l, w - 2 * p))
                             - log2(binom(k1, p) ** 2) - solutions, 0.) # permutations

                #Tp = max(log2(binom(n, w)) - log2(binom(n - k, w - 2 * p)) - log2(binom(k1, p) ** 2) - solutions, 0.)
                Tg = _gaussian_elimination_complexity2(ncols, nrows, r)
                tmp = Tp + log2(Tg + _list_merge_complexity(L1, l, hmap))

                tmp += __memory_access_cost(tmp_mem, memory_access)

                time = min(time, tmp)
                if tmp == time:
                    memory = tmp_mem
                    params = [p, l, Tp]

        for i in range(len(i_val)):
            if params[i] == i_val[i] - 1:
                stop = False
                i_val[i] += i_val_inc[i]

        if stop:
            break

    par = {"l": params[1], "p": params[0], "c": 0}
    res = {"time": time, "perms": params[2], "memory": memory, "parameters": par}
    #print(res)
    return res

# Calculated the number of expected permutations given an Instance along with the
# optimization parameters l and p.
def dumerComplextityRandomPerm(I: Instance, l: int, p: int, log = True):
    Id_size = I.nrows-l
    rs_half = math.ceil((I.ncols-Id_size)/2) #math.ceil(I.ncols - Id_size / 2)
    combs = math.ceil(  math.comb(I.ncols, I.w) / ( math.comb(Id_size, I.w - 2*p) * (math.comb(rs_half, p)**2)) )
    ret = math.log2(combs) if log else combs
    return ret

# Well, this is certainly a TODO.
def dumerComplexitySpecialPerm(I: ReducedInstance, l: int, p: int):
    pass

# Function that prints the memory consumption alongside other usefull
# information for a given ReducedInstance.
def dumerMemoryConsumption(I: ReducedInstance, l: int, p: int):
    print("Calculating and Printing Memory Consumption for params: l={}, p={}".format(l, p))
    bytes_per_element = 4
    bucket_factor = 5
    I.print_info()
    nc = I.ncols
    nr = I.nrows
    right_size = nc - (nr - l)
    right_comb_half = int(math.ceil(right_size / 2))
    nr_combs = math.comb(right_comb_half, p)
    l_table_size = 2**l
    elements_per_bucket = nr_combs / l_table_size
    bucket_size = math.ceil(elements_per_bucket) * bucket_factor
    exp_perm = dumerComplextityRandomPerm(I, l, p)
    exp_coll = log2(l_table_size *  elements_per_bucket**2)

    mem_list = bucket_size * l_table_size * bytes_per_element #in bytes
    mem_list_mb = mem_list / (1000 * 1000)
    mem_list_gb = mem_list_mb / 1000
    mem_ctr = l_table_size * 1 # just an uint8, 1 byte should be enough for the counter
    mem_ctr_mb = mem_ctr / (1000 * 1000)
    mem_ctr_gb = mem_ctr_mb / 1000
    mem_cmb_array = nr_combs * p * 2 # fix uint16
    mem_cmb_array_mb = mem_cmb_array / (1000 * 1000)
    mem_cmb_array_gb = mem_cmb_array_mb / 1000

    mem_list_mb = round(mem_list_mb,2)
    mem_list_gb = round(mem_list_gb,2)
    mem_ctr_mb = round(mem_ctr_mb,2)
    mem_ctr_gb = round(mem_ctr_gb,2)
    mem_cmb_array_mb = round(mem_cmb_array_mb,2)
    mem_cmb_array_gb = round(mem_cmb_array_gb,2)

    print("Right comb half: ", right_comb_half)
    print("Number combination: ", nr_combs)
    print("Number buckets: ", l_table_size)
    print("Elements per bucket: ", elements_per_bucket)
    print("bucket_factor: ", bucket_factor)
    print("Bucket size: ", bucket_size)
    print("mem_list: {}MB = {}GB".format(mem_list_mb, mem_list_gb))
    print("mem_ctr: {}MB = {}GB".format(mem_ctr_mb, mem_ctr_gb))
    print("mem_cmb_array: {}MB = {}GB".format(mem_cmb_array_mb, mem_cmb_array_gb))
    print("Exp. coll.: 2^", exp_coll)
    print("Exp. perm.: 2^", exp_perm)

if __name__ == '__main__':

    for k, v in MapFromNToTemplateInstance.items():
        res2 = dumer_complexity2(v)
        # l, p = res2["parameters"]["l"], res2["parameters"]["p"]
        # res = dumerComplextityRandomPerm(v, l, p)
        # resPrange = prangeComplexitySpecialPerm(v)
        print("{}, {}".format(k, res2))
        # print("{}, {}".format(k, res))
        # print("{}, {}".format(k, resPrange))
        # print(res2)

    dumerMemoryConsumption(inst2197R2, 16, 2)

    #inst156.print_info()
    #print(prangeComplexityRandomPerm(inst156))
    # print(prangeComplexityRandomPerm(inst156R))
    # print(prangeComplexityRandomPerm(inst240R))
    # print(prangeComplexityRandomPerm(inst381R))
    # print(prangeComplexityRandomPerm(inst482R))
    # print(prangeComplexityRandomPerm(inst640R))
    # print(prangeComplexityRandomPerm(inst751R))
    # print(prangeComplexityRandomPerm(inst808R))
    # print(prangeComplexityRandomPerm(inst923R))
    # print(prangeComplexityRandomPerm(inst982R))
    # print(prangeComplexityRandomPerm(inst1041R))
    # print(prangeComplexityRandomPerm(inst1101R))
    # print(prangeComplexityRandomPerm(inst1223R))

    # print("---------------------------------------------------------------------------")
    # print(prangeComplexityRandomPerm(inst156))
    #print(prangeComplexityRandomPerm(inst240))
    # print(prangeComplexityRandomPerm(inst381))
    # print(prangeComplexityRandomPerm(inst482))
    # print(prangeComplexityRandomPerm(inst640))
    # print(prangeComplexityRandomPerm(inst751))
    # print(prangeComplexityRandomPerm(inst808))
    # print(prangeComplexityRandomPerm(inst923))
    # print(prangeComplexityRandomPerm(inst982))
    # print(prangeComplexityRandomPerm(inst1041))
    # print(prangeComplexityRandomPerm(inst1101))
    # print(prangeComplexityRandomPerm(inst1223))
    #prangeComplexitySpecialPerm(inst240R)
    #prangeComplexitySpecialPerm(inst1041R)
