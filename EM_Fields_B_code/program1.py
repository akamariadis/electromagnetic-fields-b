import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

gmin, gmax_range, Ng = 0.1, 40, 400
gammas = np.linspace(gmin, gmax_range, Ng)
xs = np.logspace(-6, 2, 4000)
gmax_vals = np.full(gammas.shape, np.nan)
for i, g in enumerate(gammas):
    def g_of_x(x):
        return x * (1 + g / (1 + x**2))

    def dg_dx(x):
        return 1 + g / (1 + x**2) - (2 * g * x**2) / (1 + x**2)**2

    dvals = dg_dx(xs)
    sign_changes = np.where(dvals[:-1] * dvals[1:] < 0)[0]
    roots_crit = []

    for k in sign_changes:
        a_br, b_br = xs[k], xs[k + 1]
        try:
            rt = brentq(dg_dx, a_br, b_br)
            if rt > 0 and all(abs(rt - r) > 1e-8 for r in roots_crit):
                roots_crit.append(rt)
        except Exception:
            pass

    roots_crit.sort()
    if roots_crit:
        gvals = g_of_x(np.array(roots_crit))
        gmax_vals[i] = np.max(gvals)
mask = (gammas > 8) & ~np.isnan(gmax_vals)
G  = gammas[mask]
Yb = gmax_vals[mask] ** 2
fig, ax = plt.subplots(facecolor='white')
ax.set_facecolor('white')
ax.grid(True)
if G.size > 0:
    ax.fill_between(G, 0, Yb, color=(0.95, 0.85, 0.55), alpha=0.35)
    ax.plot(G, Yb, linewidth=2)
ymax = float(np.max(Yb)) * 1.05 if Yb.size > 0 else 1.0
ax.axvline(8, linestyle='--', linewidth=1.2)
ax.set_xlabel(r'$\gamma$')
ax.set_ylabel(r'$\kappa E_0^2$')
ax.set_title(r'Μέρος (Α): Περιοχή αστάθειας στο επίπεδο $(\gamma,\, \kappa E_0^2)$')
ax.set_xlim(gmin, gmax_range)
ax.set_ylim(0, ymax)
plt.tight_layout()
gamma = 20
kappa = 20
a     = 1

def fE(E):
    return E * (1 + gamma / (1 + kappa * E**2))

def dfdE(E):
    return (1 + gamma / (1 + kappa * E**2)
            - (2 * gamma * kappa * E**2) / (1 + kappa * E**2)**2)

Escan = np.logspace(-6, 1, 4000)
dval  = dfdE(Escan)

sign_changes = np.where(dval[:-1] * dval[1:] < 0)[0]
roots_crit = []
for k in sign_changes:
    a_br, b_br = Escan[k], Escan[k + 1]
    try:
        rt = brentq(dfdE, a_br, b_br)
        if rt > 0 and all(abs(rt - r) > 1e-8 for r in roots_crit):
            roots_crit.append(rt)
    except Exception:
        pass

roots_crit.sort()

if len(roots_crit) < 2:
    raise RuntimeError("Δεν βρέθηκαν δύο κρίσιμα σημεία. Ελέγξτε τις παραμέτρους.")

fvals       = fE(np.array(roots_crit))
f_local_min = np.min(fvals)
f_local_max = np.max(fvals)
ra = np.linspace(0.01, 0.999, 600)
E0_min = f_local_min * ra**2
E0_max = f_local_max * ra**2
E0sq_min = E0_min**2
E0sq_max = E0_max**2
fig, ax = plt.subplots(facecolor='white')
ax.set_facecolor('white')
ax.grid(True)
ax.fill_between(ra, E0sq_min, E0sq_max, color=(0.95, 0.85, 0.55), alpha=0.35)
ax.plot(ra, E0sq_min, '--k', linewidth=1.3)
ax.plot(ra, E0sq_max, '--k', linewidth=1.3)
ax.set_xlabel(r'$r/a$')
ax.set_ylabel(r'$E_0^2\ (\mathrm{V}^2/\mathrm{m}^2)$')
ax.set_title(r'Μέρος (Β): Περιοχή αστάθειας σε $(E_0^2,\ r/a)$ για $\gamma=20$, $\kappa=20$')
ax.set_xlim(0, 1)
ax.set_ylim(0, float(np.max(E0sq_max)) * 1.05)
plt.tight_layout()
gamma, kappa, a = 20, 20, 1
E0min, E0max_val, NE0 = 0.0, 3.0, 800
Emax, Ngrid = 12.0, 6000
r_over_a_list = [0.5, 0.6, 0.75, 0.95]
def fE(E):
    return E * (1 + gamma / (1 + kappa * E**2))
E0_vals = np.linspace(E0min, E0max_val, NE0)

for ra in r_over_a_list:
    S_vals   = E0_vals * (a / ra)**2
    allRoots = [None] * len(S_vals)
    Es = np.linspace(1e-12, Emax, Ngrid)
    for i, S in enumerate(S_vals):
        vals       = fE(Es) - S
        sign_changes = np.where(vals[:-1] * vals[1:] < 0)[0]
        roots_list = []
        for k in sign_changes:
            try:
                rt = brentq(lambda x, _S=S: fE(x) - _S, Es[k], Es[k + 1])
                if rt > 0:
                    roots_list.append(rt)
            except Exception:
                pass
        if roots_list:
            roots_list = np.sort(roots_list)
            keep = np.ones(len(roots_list), dtype=bool)
            for j in range(1, len(roots_list)):
                if abs(roots_list[j] - roots_list[j - 1]) < 1e-8:
                    keep[j] = False
            allRoots[i] = roots_list[keep]
        else:
            allRoots[i] = np.array([])
    inc  = np.full(len(S_vals), np.nan)

    prev = np.nan
    for i, rts in enumerate(allRoots):
        if rts is None or len(rts) == 0:
            prev = np.nan
        else:
            choice = rts[0] if np.isnan(prev) else rts[np.argmin(np.abs(rts - prev))]
            inc[i] = choice
            prev   = choice

    dec  = np.full(len(S_vals), np.nan)
    prev = np.nan
    for i in range(len(S_vals) - 1, -1, -1):
        rts = allRoots[i]
        if rts is None or len(rts) == 0:
            prev = np.nan
        else:
            choice = rts[-1] if np.isnan(prev) else rts[np.argmin(np.abs(rts - prev))]
            dec[i] = choice
            prev   = choice

    fig, ax = plt.subplots(facecolor='white')
    ax.set_facecolor('white')
    ax.grid(True)

    for i, rts in enumerate(allRoots):
        if rts is not None and len(rts) > 0:
            ax.plot(np.full_like(rts, E0_vals[i]), rts, '.k', markersize=4)

    ax.plot(E0_vals, inc, '-',  linewidth=2.0, label='ανιούσα')
    ax.plot(E0_vals, dec, '--', linewidth=2.0, label='κατιούσα')

    ax.set_xlabel(r'$E_0\ (\mathrm{V/m})$')
    ax.set_ylabel(r'$E(r)\ (\mathrm{V/m})$')
    ax.set_title(rf'Μέρος (Γ): $E(r)$ vs $E_0$ στο $r={ra:.2f}\,a$')
    ax.set_xlim(E0min, E0max_val)
    valid = np.concatenate([inc[~np.isnan(inc)], dec[~np.isnan(dec)]])
    ymax  = max(3.0, float(np.max(valid)) * 1.1) if valid.size > 0 else 3.0
    ax.set_ylim(0, ymax)
    ax.legend(loc='upper left')
    plt.tight_layout()

gamma, kappa, a = 20, 20, 1
E0_list = [0.5, 1.0, 2.0, 3.0]
rmin, rmax, Nr = 0.01, 1.5, 700
Emax, Ngrid    = 12.0, 6000

def fE(E):
    return E * (1 + gamma / (1 + kappa * E**2))
ra_vals = np.linspace(rmin, rmax, Nr)
for E0 in E0_list:
    S_vals   = E0 * (a / ra_vals)**2
    allRoots = [None] * Nr
    Es = np.linspace(1e-12, Emax, Ngrid)
    for i, S in enumerate(S_vals):
        vals         = fE(Es) - S
        sign_changes = np.where(vals[:-1] * vals[1:] < 0)[0]
        roots_list   = []
        for k in sign_changes:
            try:
                rt = brentq(lambda x, _S=S: fE(x) - _S, Es[k], Es[k + 1])
                if rt > 0:
                    roots_list.append(rt)
            except Exception:
                pass
        if roots_list:
            roots_list = np.sort(roots_list)
            keep = np.ones(len(roots_list), dtype=bool)
            for m in range(1, len(roots_list)):
                if abs(roots_list[m] - roots_list[m - 1]) < 1e-8:
                    keep[m] = False
            allRoots[i] = roots_list[keep]
        else:
            allRoots[i] = np.array([])

    inc  = np.full(Nr, np.nan)
    prev = np.nan
    for i, rts in enumerate(allRoots):
        if rts is None or len(rts) == 0:
            prev = np.nan
        else:
            choice = rts[-1] if np.isnan(prev) else rts[np.argmin(np.abs(rts - prev))]
            inc[i] = choice
            prev   = choice

    dec  = np.full(Nr, np.nan)
    prev = np.nan
    for i in range(Nr - 1, -1, -1):
        rts = allRoots[i]
        if rts is None or len(rts) == 0:
            prev = np.nan
        else:
            choice = rts[0] if np.isnan(prev) else rts[np.argmin(np.abs(rts - prev))]
            dec[i] = choice
            prev   = choice

    fig, ax = plt.subplots(facecolor='white')
    ax.set_facecolor('white')
    ax.grid(True)
    ax.plot(ra_vals, inc, '-',  linewidth=2.0, label='increasing')
    ax.plot(ra_vals, dec, '--', linewidth=2.0, label='decreasing')
    ax.set_xlabel(r'$r/a$')
    ax.set_ylabel(r'$E(r)\ (\mathrm{V/m})$')
    ax.set_title(rf'Μέρος (Δ): $E(r)$ vs $r/a$  ($E_0={E0:.1f}$ V/m)')
    ax.set_xlim(rmin, rmax)
    valid = np.concatenate([inc[~np.isnan(inc)], dec[~np.isnan(dec)]])
    ymax  = max(3.0, float(np.max(valid)) * 1.1) if valid.size > 0 else 3.0
    ax.set_ylim(0, ymax)
    ax.legend(loc='upper left')
    plt.tight_layout()

plt.show()