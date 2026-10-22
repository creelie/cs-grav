"""Numerical verification of the statements of the paper.

Run:  python3 verify.py      (about a minute on a laptop)

Every check prints PASS or FAIL; the script exits with status 1 if any
check fails.  The numerical statements are double-precision tests of
identities whose proofs are given in the paper; nothing in the paper rests
on them.
"""
import math
import sys

import numpy as np
import scipy.sparse as sp

from clifford import Clifford, antisym_product

rng = np.random.default_rng(20261022)
RESULTS = []


def report(name, ok, detail=""):
    RESULTS.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def crand(n):
    return rng.standard_normal(n) + 1j * rng.standard_normal(n)


# ---------------------------------------------------------------------------
# 1. Symmetric square of a half-spinor of so(2r): Sym^2 = sum_j Lambda^{r-4j}
# ---------------------------------------------------------------------------
def sym2_forms(r):
    return [p for p in range(r, -1, -4)]


def check_sym2():
    ok = True
    for r in range(2, 21):
        N = 2 ** (r - 1)
        lhs = N * (N + 1) // 2
        rhs = sum(math.comb(2 * r, p) // (2 if p == r else 1) for p in sym2_forms(r))
        anti = sum(math.comb(2 * r, p) for p in range(r - 2, -1, -4))
        ok &= (lhs == rhs) and (N * (N - 1) // 2 == anti)
    report("Sym^2 and Lambda^2 of the half-spinor of so(2r), r=2..20", ok)
    # vector in Sym^2 iff r = 1 mod 4, i.e. D = 2 mod 8
    ok2 = all((1 in sym2_forms(r)) == (r % 4 == 1) for r in range(2, 41))
    report("Lambda^1 in Sym^2(Delta_+) iff D = 2 mod 8", ok2)


# ---------------------------------------------------------------------------
# 2. The four families of superalgebras of Rios-Marrani-Chester are Sym^2
# ---------------------------------------------------------------------------
def check_rmc_algebras():
    ok = True
    rows = []
    for n in range(0, 4):
        fam = {
            "(9+8n)+1": (9 + 8 * n, 1, list(range(1, 6 + 4 * n, 4))),
            "(10+8n)+2": (10 + 8 * n, 2, list(range(2, 7 + 4 * n, 4))),
            "(11+8n)+3": (11 + 8 * n, 3, list(range(3, 8 + 4 * n, 4))),
            "(12+8n)+4": (12 + 8 * n, 4, list(range(0, 9 + 4 * n, 4))),
        }
        for name, (s, t, forms) in fam.items():
            r = (s + t) // 2
            ok &= sorted(sym2_forms(r)) == sorted(forms)
            rows.append((n, s, t, sorted(forms)))
    report("RHS of eqs (3.5)-(3.8) of arXiv:1811.06101 = all of Sym^2(Delta_+)", ok)
    return rows


# ---------------------------------------------------------------------------
# 3. Vector bilinear, cubic identity and the count N = 2(d-2)
# ---------------------------------------------------------------------------
def vector_blocks(c, chir_eps=+1, chir_lam=+1):
    C, sign = c.charge_conjugation(chir_eps, chir_lam)
    if C is None:
        return None
    P = {+1: c.plus, -1: c.minus}
    return [(C @ g)[P[chir_eps]][:, P[chir_lam]].tocsr() for g in c.G]


def check_cubic(s, t, expect):
    c = Clifford(s, t)
    G = vector_blocks(c, +1, +1)
    if G is None:
        report(f"cubic identity ({s},{t}): no same-chirality vector bilinear", expect is False)
        return
    sym = max(abs(g - g.T).max() if (g - g.T).nnz else 0 for g in G)
    x = crand(len(c.plus))
    Y = sum(c.eta[m] * (x @ (G[m] @ x)) * (G[m] @ x) for m in range(c.D))
    rel = np.linalg.norm(Y) / np.linalg.norm(x) ** 3
    holds = rel < 1e-10
    N = len(c.plus)
    report(f"cubic identity in ({s},{t}): holds={holds}, N={N}, 2(d-2)={2*(c.D-2)}",
           holds == expect and sym < 1e-12 and (holds == (N == 2 * (c.D - 2))),
           f"|Y|/|x|^3={rel:.2e}")


def check_R(s, t, expect):
    """Weak form R: (G_mu x)(G^mu x)^T = j G_sigma (x G^sigma x) for some j."""
    c = Clifford(s, t)
    G = vector_blocks(c, +1, +1)
    x = crand(len(c.plus))
    L = sum(c.eta[m] * np.outer(G[m] @ x, G[m] @ x) for m in range(c.D))
    Rm = sum(c.eta[m] * (x @ (G[m] @ x)) * G[m].toarray() for m in range(c.D))
    j = np.vdot(Rm.ravel(), L.ravel()) / np.vdot(Rm.ravel(), Rm.ravel())
    res = np.linalg.norm(L - j * Rm) / np.linalg.norm(L)
    holds = res < 1e-10
    report(f"weak identity R in ({s},{t}): holds={holds}" + (f", j={j.real:+.4f}" if holds else ""),
           holds == expect, f"residual {res:.2e}")


def check_traces(s, t):
    """tr(G^mu Gt_nu) = N delta and G^mu Gt_nu G_mu = -(d-2) G_nu, where Gt is
    the inverse structure Sigma* -> Sigma induced by Gamma^mu."""
    c = Clifford(s, t)
    C, _ = c.charge_conjugation(+1, +1)
    P, M = c.plus, c.minus
    Cinv = C.conj().T
    # gamma^mu : Sigma -> Sigma* is C Gamma^mu restricted; identify Sigma* with
    # Delta_- through C: then gamma^mu = Gamma^mu|_{+ -> -}, gt^mu = Gamma^mu|_{- -> +}
    g = [c.G[m][M][:, P].toarray() for m in range(c.D)]
    gt = [c.G[m][P][:, M].toarray() for m in range(c.D)]
    N = len(P)
    err1 = max(abs(np.trace(g[a] @ gt[b]) - N * c.eta[a] * (a == b)) for a in range(c.D) for b in range(c.D))
    err2 = 0.0
    for nu in range(c.D):
        S = sum(c.eta[m] * g[m] @ gt[nu] @ g[m] for m in range(c.D))
        err2 = max(err2, abs(S + (c.D - 2) * g[nu]).max())
    report(f"trace identities in ({s},{t})", err1 < 1e-9 and err2 < 1e-9, f"{err1:.1e}, {err2:.1e}")


# ---------------------------------------------------------------------------
# 4. Null-projected systems
# ---------------------------------------------------------------------------
def null_frame(c, k, kind="real"):
    """k pairs (n_i, m_i) with n_i.m_j = delta_ij, all others zero.
    kind='real': pair a timelike with a spacelike generator (needs k <= min(s,t));
    kind='complex': pair two spacelike generators with a factor i."""
    ns, ms, used = [], [], []
    sp_idx = list(range(c.t, c.D))
    tm_idx = list(range(c.t))
    for i in range(k):
        if kind == "real":
            a, b = tm_idx[i], sp_idx[i]
            ns.append((c.G[b] + c.G[a]) / math.sqrt(2))
            ms.append((c.G[b] - c.G[a]) / math.sqrt(2))
            used += [a, b]
        else:
            a, b = sp_idx[2 * i], sp_idx[2 * i + 1]
            ns.append((c.G[a] + 1j * c.G[b]) / math.sqrt(2))
            ms.append((c.G[a] - 1j * c.G[b]) / math.sqrt(2))
            used += [a, b]
    resid = [m for m in range(c.D) if m not in used]
    return ns, ms, resid


def kernel_basis(c, ns, chir):
    """Orthonormal basis of K = {lam in Delta_chir : n_i lam = 0 all i}."""
    k = len(ns)
    dimK = 2 ** (c.r - 1 - k)
    src = c.plus if (k % 2 == 0) == (chir == +1) else c.minus
    Pprod = sp.identity(c.dim, dtype=complex, format="csr")
    for n in ns:
        Pprod = Pprod @ n
    tgt = c.plus if chir == +1 else c.minus
    blk = Pprod[tgt][:, src]
    Z = blk @ (rng.standard_normal((len(src), dimK + 6)) + 0j)
    U, S, _ = np.linalg.svd(Z, full_matrices=False)
    rank = int((S > 1e-8 * S[0]).sum())
    basis = np.zeros((c.dim, rank), dtype=complex)
    basis[tgt, :] = U[:, :rank]
    # verify annihilation
    err = max(abs(n @ basis).max() for n in ns) if ns else 0.0
    return basis, rank, dimK, err


def absorbable(s, t, k, kind="real", trials=3):
    """Is f(eps gamma^mu lam)gamma_mu lam = gamma^{sigma[k]} eps v_[k] J_sigma
    solvable for J (independent of eps), with lam in K?  Returns residual."""
    c = Clifford(s, t)
    ns, ms, resid = null_frame(c, k, kind)
    lam_chir = +1
    eps_chir = +1 if c.r % 2 == 1 else -1
    C, _ = c.charge_conjugation(eps_chir, lam_chir)
    basis, rank, dimK, err = kernel_basis(c, ns, lam_chir)
    x = basis @ crand(rank)
    # full metric in the orthonormal basis: sum_mu eta_mu (eps C G^mu x) G^mu x
    Gx = [g @ x for g in c.G]
    # RHS operators Q_sigma = antisym(G_sigma, n_1, ..., n_k)
    Q = [antisym_product([c.G[sig]] + ns) if k else c.G[sig] for sig in range(c.D)]
    Pe = c.plus if eps_chir == +1 else c.minus
    A_rows, b_rows = [], []
    for _ in range(trials):
        eps = np.zeros(c.dim, dtype=complex)
        eps[Pe] = crand(len(Pe))
        lhs = sum(c.eta[m] * (eps @ (C @ Gx[m])) * Gx[m] for m in range(c.D))
        A_rows.append(np.column_stack([q @ eps for q in Q]))
        b_rows.append(lhs)
    A = np.vstack(A_rows)
    b = np.concatenate(b_rows)
    J, *_ = np.linalg.lstsq(A, b, rcond=None)
    res = np.linalg.norm(A @ J - b) / max(np.linalg.norm(b), 1e-300)
    # residual vector bilinear nonvanishing (axiom N3)
    eps = np.zeros(c.dim, dtype=complex)
    eps[Pe] = crand(len(Pe))
    dA = max(abs(eps @ (C @ Gx[m])) for m in resid) if resid else 0.0
    return res, rank, dimK, err, dA, c.D - 2 * k


def check_null(s, t, k, kind, expect):
    res, rank, dimK, err, dA, Dp = absorbable(s, t, k, kind)
    ok_dims = (rank == dimK) and err < 1e-9
    holds = res < 1e-9
    report(f"null-projected ({s},{t}), k={k} {kind} frame, residual D'={Dp}: "
           f"dim K={rank}, absorbable={holds}",
           ok_dims and holds == expect and dA > 1e-9, f"residual {res:.1e}")


# ---------------------------------------------------------------------------
# 5. Reality: a Majorana spinor annihilated by a complex null vector is zero
# ---------------------------------------------------------------------------
def majorana_B(c):
    """Find B with B conj(Gamma^mu) B^{-1} = Gamma^mu and B conj(B) = 1.
    In the Pauli basis every Gamma^mu is real or imaginary; B is the product
    of all imaginary ones (if their number is even) or of all real ones (if
    their number is odd)."""
    real = [g for g in c.G if abs(g - g.conj()).max() < 1e-12]
    imag = [g for g in c.G if abs(g + g.conj()).max() < 1e-12]
    assert len(real) + len(imag) == c.D
    for group in (imag, real):
        B = sp.identity(c.dim, dtype=complex, format="csr")
        for g in group:
            B = B @ g
        Binv = sp.csr_matrix(np.linalg.inv(B.toarray())) if c.dim <= 1024 else None
        if Binv is None:
            continue
        if all(abs(B @ g.conj() @ Binv - g).max() < 1e-9 for g in c.G):
            if abs(B @ B.conj() - sp.identity(c.dim)).max() < 1e-9:
                return B.tocsr(), Binv
    return None, None


def check_reality(s, t, k):
    c = Clifford(s, t)
    B, Binv = majorana_B(c)
    if B is None:
        report(f"Majorana structure in ({s},{t})", False)
        return
    ns, ms, _ = null_frame(c, k, "complex")
    conj_ok = all(abs(B @ n.conj() @ Binv - m).max() < 1e-9 for n, m in zip(ns, ms))
    anti = all(abs(n @ m + m @ n - 2 * sp.identity(c.dim)).max() < 1e-9 for n, m in zip(ns, ms))
    report(f"Majorana structure in ({s},{t}); conj(n_i) = m_i; n m + m n = 2",
           B is not None and conj_ok and anti)


def check_witt():
    """A real frame of k = min(s,t) mutually orthogonal null pairs exists
    (the bound k <= min(s,t) itself is Lemma 2.3 of the paper)."""
    ok = True
    for s, t in [(9, 1), (10, 2), (11, 3), (12, 4), (17, 1), (19, 3), (5, 5)]:
        c = Clifford(s, t)
        k = min(s, t)
        ns, ms, _ = null_frame(c, k, "real")
        I = sp.identity(c.dim)
        for i in range(k):
            for j in range(k):
                nn = ns[i] @ ns[j] + ns[j] @ ns[i]
                nm = ns[i] @ ms[j] + ms[j] @ ns[i]
                ok &= abs(nn).max() < 1e-12 if nn.nnz else True
                ok &= abs(nm - 2 * (i == j) * I).max() < 1e-12
    report("real null frames with k = min(s,t) pairs: n.n = 0, n_i.m_j = delta_ij", ok)


# ---------------------------------------------------------------------------
# 6. Massless multiplets: Clifford count in (9+8n,1)
# ---------------------------------------------------------------------------
def check_counting():
    ok = True
    rows = []
    for n in range(0, 4):
        D = 10 + 8 * n
        S = 2 ** (D // 2 - 1)        # real components of a MW spinor
        active = S // 2
        log2_min = active // 2       # complex dimension 2^(active/2)
        bos, fer = D - 2, S // 2
        rows.append((n, D, S, active, log2_min, bos, fer))
        if n == 0:
            ok &= (log2_min == 4 and bos == fer == 8)
        else:
            ok &= (bos != fer) and (2 ** log2_min > bos + fer)
    report("massless multiplet bound 2^{S/4} vs vector multiplet in (9+8n,1)", ok)
    return rows


def main():
    check_sym2()
    check_rmc_algebras()
    for st, ex in [((9, 1), True), ((5, 5), True), ((1, 9), True), ((17, 1), False),
                   ((13, 5), False), ((9, 9), False), ((25, 1), False)]:
        check_cubic(*st, ex)
    for st, ex in [((9, 1), True), ((5, 5), True), ((17, 1), False), ((25, 1), False)]:
        check_R(*st, ex)
    for st in [(9, 1), (17, 1)]:
        check_traces(*st)
    # Nishino's systems (s - t = 8) and the transplanted ones (s - t >= 16)
    for args, ex in [((9, 1, 0, "real"), True), ((10, 2, 1, "real"), True),
                     ((11, 3, 2, "real"), True), ((12, 4, 3, "real"), True),
                     ((13, 5, 4, "real"), True),
                     ((6, 6, 1, "real"), True), ((7, 7, 2, "real"), True),
                     ((17, 1, 0, "real"), False), ((18, 2, 1, "real"), False),
                     ((19, 3, 0, "real"), False), ((19, 3, 2, "real"), False),
                     ((20, 4, 1, "real"), False), ((20, 4, 3, "real"), False),
                     ((17, 1, 4, "complex"), True), ((19, 3, 2, "complex"), False)]:
        check_null(*args, ex)
    check_reality(17, 1, 4)
    check_witt()
    check_counting()
    nfail = RESULTS.count(False)
    print(f"\n{len(RESULTS) - nfail} checks passed, {nfail} failed")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
