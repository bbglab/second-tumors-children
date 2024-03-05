from itertools import product


cb = dict(zip('ACGT','TGCA'))


def rev(seq):
    """reverse complement of seq"""

    return ''.join(list(map(lambda s: cb[s], seq[::-1])))


def sum_dict(d1, d2):

    keys = set(list(d1.keys()) + list(d2.keys()))
    return {k: d1.get(k,0) + d2.get(k,0) for k in keys}


def triplets():

    for ref in 'CT':
        for a, b in product(cb, repeat=2):
            yield a + ref + b

def triplet_index(triplet):

    a, ref, b = tuple(list(triplet))
    s = 16 * (ref == 'T')
    t = 4 * ((a == 'C') + 2 * (a == 'G') + 3 * (a == 'T'))
    u = (b == 'C') + 2 * (b == 'G') + 3 * (b == 'T')
    return s + t + u


index_dict = {}

for i, t in enumerate(triplets()):

    assert(triplet_index(t) == i)


def sbs_format(triplet_count):
    """Maps ref triplets to 96 SBS channel"""

    vector = []
    for ref in 'CT':
        for alt in 'ACGT':
            if alt != ref:
                for a, b in product(cb, repeat=2):
                    vector.append(triplet_count[triplet_index(a + ref + b)])
    return vector


def sbs_normalize(sbs_profile, triplet_count):
    """
    sbs_profile: 96 array
    triplet_count: 32 or 96 array
    """

    if len(triplet_count) == 32:
        triplet_count = sbs_format(triplet_count)
    a = sbs_profile / triplet_count
    return a / np.sum(a)


def sbs_index(triplet, alt):

    a, ref, b = tuple(list(triplet))
    u = 48 * (ref == 'T')
    if u > 0:
        v = 16 * ((alt == 'C') + 2 * (alt == 'G'))
    else:
        v = 16 * ((alt == 'G') + 2 * (alt == 'T'))
    s = 4 * ((a == 'C') + 2 * (a == 'G') + 3 * (a == 'T'))
    t = (b == 'C') + 2 * (b == 'G') + 3 * (b == 'T')
    return u + v + s + t

def count_triplets(seq):

    return [seq.count(t) + seq.count(rev(t)) for t in triplets()]


def mut_key_gen():

    for ref in ['C', 'T']:
        for alt in cb.keys():
            if ref == alt:
                continue
            else:
                for p in product(cb.keys(), repeat=2):
                    yield p[0] + ref + p[1], alt

# if __name__ == '__main__':
# 
    # for triplet, alt in mut_key_gen():
        # print(triplet +'>'+alt)
        # print(sbs_index(triplet, alt))
