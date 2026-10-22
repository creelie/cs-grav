"""Explicit Clifford algebras Cl(s,t) with s+t even, built from Pauli matrices.

All matrices are scipy sparse.  The first t generators are timelike
(square -1), the remaining s are spacelike (square +1).  The chirality
operator is diagonal in this basis, so the Weyl subspaces are coordinate
subspaces.
"""
import itertools
import math

import numpy as np
import scipy.sparse as sp

S1 = sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
S2 = sp.csr_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex))
S3 = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))
I2 = sp.identity(2, dtype=complex, format="csr")


def kron_list(mats):
    out = mats[0]
    for m in mats[1:]:
        out = sp.kron(out, m, format="csr")
    return out


def euclidean_gammas(r):
    """2r Hermitian generators of Cl(2r,0) of size 2^r (Jordan-Wigner)."""
    gs = []
    for j in range(r):
        left = [S3] * j
        right = [I2] * (r - j - 1)
        gs.append(kron_list(left + [S1] + right) if r > 1 else S1)
        gs.append(kron_list(left + [S2] + right) if r > 1 else S2)
    return gs


class Clifford:
    def __init__(self, s, t):
        D = s + t
        assert D % 2 == 0 and D >= 2
        self.s, self.t, self.D, self.r = s, t, D, D // 2
        e = euclidean_gammas(self.r)
        self.eta = np.array([-1.0] * t + [1.0] * s)
        self.G = [(1j * e[a] if a < t else e[a]) for a in range(D)]
        self.dim = 2 ** self.r
        chi = np.ones(1)
        for _ in range(self.r):
            chi = np.kron(chi, np.array([1.0, -1.0]))
        self.chi = chi
        self.plus = np.where(chi > 0)[0]
        self.minus = np.where(chi < 0)[0]
        # charge conjugation candidates: products of the antisymmetric
        # (sigma_2-type) and of the symmetric (sigma_1-type) Euclidean gammas
        self.C_candidates = []
        for parity in (0, 1):
            C = sp.identity(self.dim, dtype=complex, format="csr")
            for a in range(self.D):
                if a % 2 == parity:
                    C = C @ e[a]
            self.C_candidates.append(C.tocsr())

    def check_clifford(self):
        err = 0.0
        for a in range(self.D):
            for b in range(self.D):
                M = self.G[a] @ self.G[b] + self.G[b] @ self.G[a]
                target = 2 * self.eta[a] * (a == b) * sp.identity(self.dim)
                err = max(err, abs(M - target).max() if (M - target).nnz else 0.0)
        return err

    def charge_conjugation(self, chir_eps, chir_lam):
        """Return C such that B(x,y)=x^T C y is Spin-invariant and the vector
        bilinear x^T C Gamma^mu y pairs chirality chir_eps with chir_lam
        nontrivially.  Returns (C, sign) with C Gamma C^{-1} = sign Gamma^T,
        or (None, None)."""
        P = {+1: self.plus, -1: self.minus}
        for C in self.C_candidates:
            Cinv = C.conj().T.tocsr()  # C is unitary
            for sign in (+1, -1):
                ok = all(abs(C @ g @ Cinv - sign * g.T).max() < 1e-9 for g in self.G)
                if not ok:
                    continue
                blk = (C @ self.G[self.t if self.s else 0])[P[chir_eps]][:, P[chir_lam]]
                if abs(blk).max() > 1e-9:
                    return C, sign
        return None, None


def antisym_product(mats):
    """Totally antisymmetrised product (1/n!) sum sign(p) M_p1 ... M_pn."""
    n = len(mats)
    out = None
    for perm in itertools.permutations(range(n)):
        sgn = np.linalg.det(np.eye(n)[list(perm)])
        prod = mats[perm[0]]
        for i in perm[1:]:
            prod = prod @ mats[i]
        out = sgn * prod if out is None else out + sgn * prod
    return out / math.factorial(n)
