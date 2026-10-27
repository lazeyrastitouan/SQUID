import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plot_upper_limit(ax, x, upper, color, label, fade_dex=1.0, alpha0=0.5, nlayers=30):
    """Solid line at the upper limit, with a shaded region below it that fades to
    transparent over `fade_dex` decades (so the y-range is not stretched to tiny values)."""
    upper = np.asarray(upper, dtype=float)
    ax.plot(x, upper, color=color, lw=1.5, label=label)
    for k in range(nlayers):
        lo = upper * 10**(-fade_dex*k/nlayers)
        hi = upper * 10**(-fade_dex*(k+1)/nlayers)
        a  = alpha0 * (1 - k/nlayers)      # linear fade from alpha0 to ~0
        ax.fill_between(x, hi, lo, color=color, alpha=a, linewidth=0)

