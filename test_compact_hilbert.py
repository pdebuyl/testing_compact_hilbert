from __future__ import print_function, division
import numpy as np
import math
import h5py
import sys
import os.path
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('M', type=int, nargs='+')
args = parser.parse_args()

filename = 'hilbert_'+'x'.join(map(str, args.M))+'.h5'

if not os.path.isfile(filename):
    print('The reference file', filename, 'does not exist')
    sys.exit(0)

with h5py.File(filename, 'r') as a:
    file_data = a['value'][:]
    assert np.all(np.equal(a['value'].attrs['M'], args.M))

N = len(args.M)
M = max(args.M)

def bin_str(i):
    """Return a string representation of i with N bits."""
    out = ''
    for j in range(N-1,-1,-1):
        if (i>>j) & 1 == 1:
            out += '1'
        else:
            out += '0'
    return out

def rotate_right(x, d):
    """Rotate x by d bits to the right."""
    d = d % N
    out = x >> d
    for i in range(d):
        bit = (x & 2**i)>>i
        out |= bit << (N+i-d)
    return out

def rotate_left(x, d):
    """Rotate x by d bits to the left."""
    d = d % N
    out = x << d
    excess = out 
    out = out & (2**N-1)
    for i in range(d):
        bit = (x & 2**(N-1-d+1+i))>> (N-1-d+1+i)
        out |= bit << i
    return out

def bit_component(x, i):
    """Return i-th bit of x"""
    return (x & 2**i) >> i


# verify that '~i & 2**N-1' performs the NOT operation in N-bit space
for i in range(2**N):
    not_i = ~i & 2**N-1
    assert not_i >=0
    assert not_i < 2**N
    assert i & not_i == 0
    assert i | not_i == 2**N-1


def gc(i):
    """Return the Gray code index of i."""
    return i ^ (i >> 1)

def e(i):
    """Return the entry point of hypercube i."""
    if i==0:
        return 0
    else:
        return gc(2*int(math.floor((i-1)//2)))

def f(i):
    """Return the exit point of hypercube i."""
    return e(2**N-1-i) ^ 2**(N-1)

def i_to_p(i):
    """Extract the 3d position from a 3-bit integer."""
    return [bit_component(i,j) for j in (0,1,2)]

def inverse_gc(g):
    """The inverse gray code."""
    i = g
    j = 1
    while j<N:
        i = i ^ (g >> j)
        j = j + 1
    return i

def g(i):
    """The direction between subcube i and the next one"""
    return int(np.log2(gc(i)^gc(i+1)))


def d(i):
    """The direction of the arrow whithin a subcube."""
    if i==0:
        return 0
    elif (i%2)==0:
        return g(i-1) % N
    else:
        return g(i) % N

def T(e, d, b):
    """Transform b."""
    out = b ^ e
    return rotate_right(out, d+1)

def T_inv(e, d, b):
    """Inverse transform b."""
    return T(rotate_right(e, d+1), N-d-2, b)


# only g(i)-th bit changes from gc(i) to gc(i+1)
for i in range(2**N-1):
    assert gc(i) ^ (1 << g(i)) == gc(i+1)

# T sends e(i) to 0 and f(i) to 2**(N-1)
for i in range(2**N):
    assert T(e(i), d(i), e(i))==0
    assert T(e(i), d(i), f(i))==2**(N-1)

# e(i) reflected in direction d(i) is f(i)
for i in range(2**N):
    assert e(i) ^ (1<<d(i)) == f(i)

# T_inv composed with T (and vice versa) is the identity operator
for i in range(2**N):
    for b in range(2**N):
        assert T(e(i), d(i), T_inv(e(i), d(i), b)) == b
        assert T_inv(e(i), d(i), T(e(i), d(i), b)) == b

def TR_algo2(p, vd=2):
    """Return the Hilbert index of point p"""
    # h will contain the Hilbert index
    h = 0
    # ve and vd contain the entry point and dimension of the current subcube
    # we choose here a main traversal direction N-2 (i.e. z for a cube) to match
    # the illustrations
    ve = 0
    for i in range(M-1, -1, -1):
        # the cell label is constructed in two steps
        # 1. extract the relevant bits from p
        l = [bit_component(px, i) for px in p]
        # 2. construct a integer whose bits are given by l
        l = sum( [lx*2**j for j, lx in enumerate(l)] )
        # transform l into the current subcube
        l = T(ve, vd, l)
        # obtain the gray code ordering from the label l
        w = inverse_gc(l)
        # compose (see [TR] lemma 2.13) the transform of ve and vd
        # with the data of the subcube
        ve = ve ^ (rotate_left(e(w), vd+1))
        vd = (vd + d(w) + 1) % N
        # move the index to more significant bits and add current value
        h = (h << N) | w
    return h

def TR_algo3(h, vd=2):
    """Return the coordinates for the Hilbert index h"""
    ve = 0
    vd = 2
    p = [0]*N
    for i in range(M-1, -1, -1):
        w = [bit_component(h, i*N+ii) for ii in range(N)]
        #print(i, w)
        w = sum( [wx*2**j for j, wx in enumerate(w)] )
        #print(i, w, gc(w))
        l = gc(w)
        l = T_inv(ve, vd, l)
        for j in range(N):
            p[j] += bit_component(l, j) << i
        ve = ve ^ rotate_left(e(w), vd+1)
        vd = (vd + d(w) + 1) % N
    return p


# Verifying that the algorithm and its inverse agree
for h_idx in range(2**(N*M)):
    assert TR_algo2(TR_algo3(h_idx)) == h_idx

def gcr(i, mu, pi):
    """Compute the gray code rank of i given the mask mu.
    Algorithm 4 in [TR]"""
    r = 0
    for k in range(N-1, -1, -1):
        if bit_component(mu, k):
            r = (r << 1) | bit_component(i, k)
    return r

def gcr_inv(r, mu, pi):
    """Inverse of the gray code rank, given the mask mu and the pattern pi.
    Algorithm 5 in [TR]"""
    i = 0
    g = 0
    j = sum([bit_component(mu, k) for k in range(N)])-1
    for k in range(N-1, -1, -1):
        if bit_component(mu, k)==1:
            i |= bit_component(r, j) << k
            g |= ( (bit_component(i, k) + bit_component(i, k+1))%2 ) << k
            j -= 1
        else:
            g |= bit_component(pi, k) << k
            i |= ( (bit_component(g, k) + bit_component(i, k+1)) % 2) << k
    return i

class CompactHilbert(object):

    def __init__(self, compact_M, vd=2):
        self._compact_M = compact_M
        self._N = len(compact_M)
        self._vd = vd
    
    @property
    def hmax(self):
        return 2**sum(self._compact_M)
    
    def extract_mask(self, i):
        """Extract the mask for iteration i of the algorithm.
        Algorithm 6 in [TR]"""
        mu = 0
        for j in range(self._N-1, -1, -1):
            mu = mu << 1
            if self._compact_M[j] > i:
                mu = mu | 1
        return mu

    def TR_algo7(self, p):
        """Compute the compact Hilbert index for point p"""
        h = 0
        ve = 0
        vd = self._vd
        m = max(self._compact_M)
        for i in range(m-1, -1, -1):
            mu = self.extract_mask(i)
            mu_norm = sum([bit_component(mu, j) for j in range(self._N)])
            mu = rotate_right(mu, vd+1)
            pi = rotate_right(ve, vd+1) & ((~mu) & 2**self._N-1)
            l = [bit_component(px, i) for px in p]
            # 2. construct a integer whose bits are given by l
            l = sum( [lx*2**j for j, lx in enumerate(l)] )
            l = T(ve, vd, l)
            w = inverse_gc(l)
            r = gcr(w, mu, pi)
            ve = ve ^ rotate_left(e(w), vd+1)
            vd = (vd + d(w) + 1) % self._N
            h = (h << mu_norm) | r
        return h

    def TR_algo8(self, h):
        """Compute the point with compact Hilbert index h"""
        ve = 0
        vd = self._vd
        k = 0
        p = [0,]*self._N
        m = max(self._compact_M)
        vM = sum(self._compact_M)
        for i in range(m-1, -1, -1):
            mu = self.extract_mask(i)
            mu_norm = sum([bit_component(mu, j) for j in range(self._N)])
            mu = rotate_right(mu, vd+1)
            pi = rotate_right(ve, vd+1) & (~mu & 2**self._N-1)
            r = [bit_component(h, vM - k - (j+1)) for j in range(mu_norm)][::-1]
            r = sum( [rx*2**j for j, rx in enumerate(r)] )
            k = k + mu_norm
            w = gcr_inv(r, mu, pi)
            l = gc(w)
            l = T_inv(ve, vd, l)
            for j in range(self._N):
                p[j] |= bit_component(l, j) << i
            ve = ve ^ (rotate_left(e(w), vd+1))
            vd = (vd + d(w) + 1) % self._N
        return p

tester = CompactHilbert(compact_M=args.M, vd=0)

test_data = [tester.TR_algo8(idx) for idx in range(tester.hmax)]
    
assert np.all(np.equal(test_data, file_data))

hilbert_cube_indices = np.array([TR_algo2(p, vd=0) for p in test_data])

assert (hilbert_cube_indices[1:]-hilbert_cube_indices[:-1]).min()>0
