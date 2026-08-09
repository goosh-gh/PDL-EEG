import numpy as np
np.random.seed(7)

# --- tiny synthetic leadfield, Ne != Ns so any transpose bug dies immediately ---
Ne, Ns = 8, 20
K = np.random.randn(Ne, Ns)          # math shape: electrodes x sources (b = K j)

def centering(n):                    # average-reference operator H = I - 11^T/n
    return np.eye(n) - np.ones((n,n))/n

def gram(K, alpha, ref):
    Ne = K.shape[0]
    H = centering(Ne) if ref=='car' else np.eye(Ne)
    C = K @ K.T + alpha*H
    return C, H

def mne_operator(K, alpha=1e-2, ref='none'):
    C,_ = gram(K, alpha, ref)
    return K.T @ np.linalg.inv(C)               # Ns x Ne

def sloreta_diag(K, alpha=1e-2, ref='none'):
    # standardizer S_ii = k_i^T C^{-1} k_i  (resolution diag), scalar sources
    C,_ = gram(K, alpha, ref)
    Ci = np.linalg.inv(C)
    T  = K.T @ Ci                                # MNE operator
    R  = T @ K                                   # full resolution (small here)
    return T, np.diag(R).copy()

def eloreta_weights(K, alpha=1e-2, ref='none', n_iter=200, tol=1e-12):
    Ne, Ns = K.shape
    H = centering(Ne) if ref=='car' else np.eye(Ne)
    w = np.ones(Ns)                              # diagonal weights (scalar sources)
    for it in range(n_iter):
        Winv = 1.0/w
        M = np.linalg.inv((K*Winv) @ K.T + alpha*H)   # (K W^{-1} K^T + aH)^{-1}
        # w_i <- sqrt(k_i^T M k_i)
        wn = np.sqrt(np.einsum('ji,jk,ki->i', K, M, K))
        rel = np.max(np.abs(wn-w))/np.max(wn)
        w = wn
        if rel < tol:
            break
    Winv = 1.0/w
    M = np.linalg.inv((K*Winv) @ K.T + alpha*H)
    T = (K.T * Winv[:,None]) @ M                 # Ns x Ne  eLORETA operator
    return T, w, it+1, rel

# ---- localization tests (noiseless, small alpha) ----
def loc_err(T, standardize=None, active=(5,), amp=None):
    j = np.zeros(Ns);  amp = amp or [1.0]*len(active)
    for a,v in zip(active, amp): j[a]=v
    b = K @ j
    jhat = T @ b
    power = jhat**2
    if standardize is not None:
        power = jhat**2 / standardize
    est = set(np.argsort(power)[-len(active):])
    return est, set(active), power

al = 1e-6
T_s, Sdiag = sloreta_diag(K, alpha=al, ref='none')
T_e, w, nit, rel = eloreta_weights(K, alpha=al, ref='none')

print("=== sLORETA single-source localization (each source) ===")
ok=0
for src in range(Ns):
    est,tru,_ = loc_err(T_s, standardize=Sdiag, active=(src,))
    ok += (est==tru)
print(f"  exact-hit {ok}/{Ns}")

print("=== eLORETA single-source ===")
ok=0
for src in range(Ns):
    est,tru,_ = loc_err(T_e, standardize=None, active=(src,))
    ok += (est==tru)
print(f"  exact-hit {ok}/{Ns}  (iters={nit}, rel={rel:.1e})")

print("=== eLORETA two-source (exact-zero-error property) ===")
ok=0; N2=0
for a in range(0,Ns,3):
    for b_ in range(a+2,Ns,5):
        est,tru,_ = loc_err(T_e, standardize=None, active=(a,b_), amp=[1.0,0.7])
        ok += (est==tru); N2+=1
print(f"  exact-hit {ok}/{N2}")

# dump matrices + a canonical reconstruction for PDL cross-check
np.save('K.npy', K)
j = np.zeros(Ns); j[5]=1.0; j[11]=-0.6
b = K@j
np.savez('canon.npz', K=K, b=b,
         Tmne=mne_operator(K,alpha=al,ref='none'),
         Ts=T_s, Sdiag=Sdiag, Te=T_e, w=w,
         jhat_mne=(mne_operator(K,alpha=al,ref='none')@b),
         jhat_slor=(T_s@b), jhat_elor=(T_e@b))
print("\nsaved canon.npz  (alpha=%.0e, iters=%d)"%(al,nit))
print("w[:5]=",np.round(w[:5],6))
