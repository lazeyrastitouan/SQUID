"""R_tot * Pfav(R) as an interpolating function.

    from rate_pfav import rate_pfav, rate_pfav_slice
    f = rate_pfav(sigma=0.6, m_min=1, network='ET_triangle', channel='PBH', z=0.3)
    f(60.)                 # Gpc^-3 yr^-1, R in comoving Mpc/h
    f.no_localisation      # value with no localisation requirement
    g = rate_pfav_slice(sigma=0.6, m_min=1, network='ET_triangle', channel='PBH')   # averaged over 0.15 < z < 0.43
"""
from pathlib import Path
import numpy as np

_T = np.genfromtxt(Path(__file__).resolve().parent / 'rates_x_pfav.csv', delimiter=',',
                   names=True, dtype=None, encoding='utf-8')
Z_GRID = np.unique(_T['z'])
KAPPA = {'PBH': 1.28, 'BBH': 2.9, 'BNS': 0.0}          # R_tot(z) proportional to (1+z)^KAPPA
OM = 0.29                                             # flat LCDM, h = 0.6831


def _rows(z, sigma, m_min, network, channel, rate_case, m_c):
    sel = ((_T['z'] == z) & (_T['sigma'] == sigma) & (_T['m_c'] == m_c) & (_T['m_min'] == m_min)
           & (_T['network'] == network) & (_T['channel'] == channel) & (_T['rate_case'] == rate_case))
    if not sel.any():
        raise ValueError(f'no entry for z={z}, sigma={sigma}, m_c={m_c}, m_min={m_min}, {network}, {channel}, {rate_case}')
    return _T[sel]


def _make(R, y, M_star):
    fin = np.isfinite(R)
    Rg, yg = R[fin], y[fin]

    def f(r):
        # linear interpolation on the 1 Mpc/h grid; 0 below 10 Mpc/h, held constant above 200 Mpc/h
        return np.interp(r, Rg, yg, left=0.0, right=yg[-1])

    f.no_localisation = float(y[~fin][0])
    f.R_grid, f.values, f.M_star = Rg, yg, M_star
    return f


def _pfav_z(z, sigma, m_min, network, channel, rate_case, m_c):
    """R grid, Pfav(R, z) (linear in z between Fisher redshifts), R_tot(z), M_star."""
    zs = Z_GRID if m_c == 1.0 else np.array([0.3])
    if not zs[0] <= z <= zs[-1]:
        raise ValueError(f'z={z} outside the Fisher redshifts for m_c={m_c}: {zs[0]:g} to {zs[-1]:g}')
    i = min(np.searchsorted(zs, z, side='right') - 1, len(zs) - 2) if len(zs) > 1 else 0
    lo = _rows(zs[i], sigma, m_min, network, channel, rate_case, m_c)
    if len(zs) == 1:
        return lo['R_mpch'], lo['Pfav'], lo['R_tot'][0], float(lo['M_star'][0])
    hi = _rows(zs[i + 1], sigma, m_min, network, channel, rate_case, m_c)
    t = (z - zs[i]) / (zs[i + 1] - zs[i])
    rtot = lo['R_tot'][0] * ((1 + z) / (1 + zs[i]))**KAPPA[channel]
    return lo['R_mpch'], (1 - t) * lo['Pfav'] + t * hi['Pfav'], rtot, float(lo['M_star'][0])


def _case(channel, rate_case):
    if channel == 'PBH':
        return 'nominal'
    if rate_case is None:
        raise ValueError(f"rate_case required for {channel}: {sorted(str(c) for c in set(_T['rate_case'][_T['channel'] == channel]))}")
    return rate_case


def rate_pfav(sigma, m_min, network, channel, rate_case=None, m_c=1.0, z=0.3):
    """channel: 'PBH', 'BBH' or 'BNS'.  rate_case: required for BBH ('lo', 'median', 'hi') and BNS ('lo', 'hi').
    z: source redshift, 0.15 to 0.43 (m_c = 1) or 0.3 only (m_c = 5)."""
    R, P, rtot, M_star = _pfav_z(z, sigma, m_min, network, channel, _case(channel, rate_case), m_c)
    return _make(R, rtot * P, M_star)


def rate_pfav_slice(sigma, m_min, network, channel, rate_case=None, m_c=1.0, z_min=0.15, z_max=0.43, nz=57):
    """Comoving-volume average over z_min < z < z_max of R_tot(z) Pfav(R, z) / (1+z)  [Gpc^-3 yr^-1, observer time].
    Number of events: V_survey [Gpc^3] x T_obs [yr] x (this, with f_v V_frac C applied as for rate_pfav)."""
    case = _case(channel, rate_case)
    zz = np.linspace(z_min, z_max, nz)
    zf = np.linspace(0, z_max, 4001)
    dc = np.concatenate([[0], np.cumsum(np.diff(zf) * 0.5 * (1 / np.sqrt(OM * (1 + zf[1:])**3 + 1 - OM)
                                                             + 1 / np.sqrt(OM * (1 + zf[:-1])**3 + 1 - OM)))])
    w = np.interp(zz, zf, dc)**2 / np.sqrt(OM * (1 + zz)**3 + 1 - OM)          # dV/dz up to a constant
    w = w * np.gradient(zz)
    w[[0, -1]] *= 0.5                                                       # trapezoid
    tot = 0
    for zi, wi in zip(zz, w):
        R, P, rtot, M_star = _pfav_z(zi, sigma, m_min, network, channel, case, m_c)
        tot = tot + wi * rtot * P / (1 + zi)
    return _make(R, tot / w.sum(), M_star)
