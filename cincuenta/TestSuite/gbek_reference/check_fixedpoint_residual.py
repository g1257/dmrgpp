#!/usr/bin/env python3
"""
One-off check (2026-07-09): compare err[ch]/err[ev] at rank L=3 for two
different candidate "targets":

  (A) cincuenta's own C++ atomic-limit reference dump
      (atomic-limit-gbek-L3-weiss-delta-lesser) -- built from the FULL
      extended-Fock/exact dynamics, NOT compressed through a rank-3 Cholesky
      loop anywhere upstream.
  (B) gbek_selfconsistency.py's own reference (gbek-atomic-limit-exact-lesser)
      -- built by ITERATING the rank-3 causal Cholesky decomposition inside
      the DMFT self-consistency loop itself (see run_self_consistency() in
      gbek_selfconsistency.py: `V = cholesky_causal(Lambda, L)` every
      iteration), i.e. a genuine fixed point of a rank-3-compressed loop,
      analogous to what GBEK's own Fig. 3 target actually is (per the
      "strongest Fable" hypothesis under test: their target is co-generated
      with rank-3 causal compression, not a rank-independent exact target).

If that hypothesis is right, target (B)'s err[ch]/err[ev] ratio should sit
much closer to the paper's own quoted 0.17/0.09 = 1.9x than target (A)'s.
"""
import numpy as np

from compare_reference import read_lesser_file
from gbek_cholesky import cholesky_causal, reconstruct

TARGET_A = "/Users/Shared/ornldev/code/dmrgpp/build/atomic-limit-gbek-L3-weiss-delta-lesser"
TARGET_B = "/tmp/reconv_check.txt"  # freshly regenerated, post seeding-floor-fix
L = 3


def load_lambda_from_delta(path):
    ts, re, im = read_lesser_file(path)
    return ts, -1j * (re + 1j * im)


def load_lambda_direct(path):
    ts, re, im = read_lesser_file(path)
    return ts, re + 1j * im


def eigenvector_lowrank(lam, L):
    w, U = np.linalg.eigh(lam)
    order = np.argsort(w)[::-1]
    w, U = w[order], U[:, order]
    return (U[:, :L] * w[:L]) @ U[:, :L].conj().T


def global_err(lam, A):
    return np.linalg.norm(lam - A) / np.linalg.norm(lam)


def report(name, ts, lam):
    V = cholesky_causal(lam, L)
    ch = reconstruct(V)
    ev = eigenvector_lowrank(lam, L)
    e_ch = global_err(lam, ch)
    e_ev = global_err(lam, ev)
    print(f"{name}: err[ch]={e_ch:.4f}  err[ev]={e_ev:.4f}  ratio={e_ch/e_ev:.2f}x  "
          f"diag(t=4)={lam[-1,-1].real:.4f}")
    return e_ch, e_ev


if __name__ == "__main__":
    tsA, lamA = load_lambda_from_delta(TARGET_A)
    tsB, lamB = load_lambda_direct(TARGET_B)

    print(f"paper (quoted, Fig. 3 L=3): err[ch]=0.17  err[ev]=0.09  ratio=1.9x\n")
    report("(A) cincuenta C++ exact (no rank-3 in loop)", tsA, lamA)
    report("(B) gbek_selfconsistency.py (rank-3 IN the loop)", tsB, lamB)
