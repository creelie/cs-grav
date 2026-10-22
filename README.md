# Supersymmetric Yang–Mills theories with several times reduce to ten dimensions

Paper, LaTeX source and verification code. Author: Deep Bhattacharjee ([ORCID 0000-0003-0466-750X](https://orcid.org/0000-0003-0466-750X)).

**Paper:** [paper/paper1/sym-several-times.pdf](https://github.com/creelie/cs-grav/blob/main/paper/paper1/sym-several-times.pdf) (JHEP format, 25 pages) · [one-column version](https://github.com/creelie/cs-grav/blob/main/paper/paper1/sym-several-times-onecolumn.pdf) (22 pages) · [LaTeX source](https://github.com/creelie/cs-grav/blob/main/paper/paper2/sym-several-times-tex.zip)

## What the paper proves

Supersymmetric Yang–Mills theories in signatures with several times are built by restricting the gaugino to the kernel of k mutually orthogonal null vectors n<sub>1</sub>, …, n<sub>k</sub> and letting the supersymmetry close on translations contracted with n<sub>1</sub> ∧ … ∧ n<sub>k</sub>. Nishino constructed such theories in the signatures (9+m)+(1+m). It has been claimed that his existence proofs carry over to (9+8n)+1, (10+8n)+2, (11+8n)+3 and (12+8n)+4, which include the exceptional signatures 27+3 and 19+3.

**Main Theorem.** Let a supersymmetric Yang–Mills theory with a Majorana–Weyl gaugino in signature (s,t) be built on a frame of k null vectors in the manner of Nishino (axioms N1–N5 of the paper). Then

- k ≤ min(s,t);
- the residual dimension s+t−2k is 10 or 2;
- s−t ∈ {8, 0, −8}.

In particular there is no such theory in (9+8n)+1, (10+8n)+2, (11+8n)+3 or (12+8n)+4 for n ≥ 1, and none in 27+3 or 19+3.

The proof reduces the cubic fermion terms to a Fierz identity on the residual spinors and shows (Theorem 3.3) that symmetric gamma matrices satisfying that identity on a d-dimensional space have spinor dimension 2(d−2), which happens only for d = 3, 4, 6, 10, the division-algebra dimensions.

**Single time (Theorem 6.1).** In (9+8n)+1 with n ≥ 1, every massless supermultiplet of a unitary theory without tensor charges has at least 2<sup>2<sup>2+4n</sup></sup> states, half of them bosonic. So a unitary theory with no tensor charges on its massless states has no massless vector multiplet: for n = 1 a vector multiplet would have 16 bosonic states, while the bound is 2<sup>64</sup>.

## Layout

| Path | Contents |
|---|---|
| `paper/paper1/` | the compiled paper: JHEP format and one-column format |
| `paper/paper2/` | `sym-several-times-tex.zip`, the LaTeX source with figures |
| `paper/tex/` | the same LaTeX source, unzipped |
| `code/clifford.py` | explicit real and complex gamma matrices in any signature |
| `code/verify.py` | 34 numerical checks of the statements of the paper |
| `code/verify_output.txt` | output of the last run |
| `code/figures/make_figures.py` | draws the plotted figures into `paper/tex/figures/` |

## Reproducing

```
pip install numpy scipy matplotlib
cd code && python3 verify.py                 # ends with "34 checks passed, 0 failed"
cd figures && python3 make_figures.py        # redraws the plotted figures
cd ../../paper/tex && latexmk -pdf main.tex   # JHEP version; hexagon.tex for one column
```

The checks are double-precision tests in explicit gamma matrices of identities that the paper proves by hand: the decomposition of Sym²Δ<sub>+</sub>, the cubic identity and the count N = 2(d−2), absorption of the cubic terms for null frames in Nishino's signatures and its failure in the families above (including 19+3), real null frames with k = min(s,t), and the state count of Theorem 6.1. No proof in the paper depends on them.

## Citing

Cite the paper and this archive. Release v1.0.0 is archived on Zenodo under [doi:10.5281/zenodo.23558751](https://doi.org/10.5281/zenodo.23558751), and each later release under a DOI of its own; the DOI [10.5281/zenodo.23558750](https://doi.org/10.5281/zenodo.23558750) stands for all versions and resolves to the latest. See `CITATION.cff`.

## License

MIT, see `LICENSE`.
