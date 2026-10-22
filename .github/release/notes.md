# v1.0.1: Supersymmetric Yang–Mills theories with several times reduce to ten dimensions

Patch release of the paper by Deep Bhattacharjee and its verification code. The mathematics and the code are those of v1.0.0 (doi:10.5281/zenodo.23558751); the paper now cites that archive by its DOI, and Table 2 is set at the end of Appendix C, before the references, instead of between references [1] and [2].

**Result.** A supersymmetric Yang–Mills theory with a Majorana–Weyl gaugino in signature (s,t), built on k mutually orthogonal null vectors in the manner of Nishino's theories in (9+m)+(1+m) (axioms N1–N5 of the paper), has k ≤ min(s,t), residual dimension s+t−2k ∈ {2,10}, and therefore s−t ∈ {8,0,−8}. There is no such theory in (9+8n)+1, (10+8n)+2, (11+8n)+3 or (12+8n)+4 for n ≥ 1, including 27+3 and 19+3. In (9+8n)+1 with n ≥ 1 every massless supermultiplet without tensor charges has at least 2^(2^(2+4n)) states, so a unitary theory with no tensor charges on its massless states has no massless vector multiplet.

**Assets.**
- `sym-several-times.pdf`: the paper in JHEP format (25 pages).
- `sym-several-times-onecolumn.pdf`: the same paper in a one-column article layout (22 pages).
- `sym-several-times-tex.zip`: LaTeX source and figures; `pdflatex main.tex` (JHEP) or `pdflatex hexagon.tex` (one column), each run three times.

The code in `code/` checks the statements in explicit gamma matrices (`python3 code/verify.py`, 34 checks) and draws the plotted figures (`python3 code/figures/make_figures.py`). Every statement of the paper is proved in the text; no proof depends on the code.
