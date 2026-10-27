"""R_tot * Pfav(R) as an interpolating function.

    from rate_pfav import rate_pfav
    f = rate_pfav(sigma=0.6, m_min=1, network='ET_triangle', channel='PBH')
    f(60.)                 # Gpc^-3 yr^-1, R in comoving Mpc/h
    f.no_localisation      # value with no localisation requirement
"""
from pathlib import Path
import numpy as np

_T = np.genfromtxt(Path(__file__).resolve().parent / 'rates_x_pfav.csv', delimiter=',',
                   names=True, dtype=None, encoding='utf-8')


def rate_pfav(sigma, m_min, network, channel, rate_case=None, m_c=1.0):
    """channel: 'PBH', 'BBH' or 'BNS'.  rate_case: required for BBH ('lo', 'median', 'hi') and BNS ('lo', 'hi')."""
    if channel == 'PBH':
        rate_case = 'nominal'
    elif rate_case is None:
        raise ValueError(f"rate_case required for {channel}: {sorted(str(c) for c in set(_T['rate_case'][_T['channel'] == channel]))}")
    sel = ((_T['sigma'] == sigma) & (_T['m_c'] == m_c) & (_T['m_min'] == m_min) & (_T['network'] == network)
           & (_T['channel'] == channel) & (_T['rate_case'] == rate_case))
    if not sel.any():
        raise ValueError(f'no entry for sigma={sigma}, m_c={m_c}, m_min={m_min}, {network}, {channel}, {rate_case}')
    R, y = _T['R_mpch'][sel], _T['R_tot_x_Pfav'][sel]
    fin = np.isfinite(R)
    Rg, yg = R[fin], y[fin]

    def f(r):
        # linear interpolation on the 1 Mpc/h grid; 0 below 10 Mpc/h, held constant above 200 Mpc/h
        return np.interp(r, Rg, yg, left=0.0, right=yg[-1])

    f.no_localisation = float(y[~fin][0])
    f.R_grid, f.values, f.M_star = Rg, yg, float(_T['M_star'][sel][0])
    return f
