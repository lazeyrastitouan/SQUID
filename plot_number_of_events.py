import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from plot_fancy import plot_upper_limit, add_y_jump   # plot_fancy.py is left untouched

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

#------------------------------------------------------------------ settings
config  = 'ET_2L'
years   = 10
yfactor = 19.53/3.61*8.02/23.22        # BOSS-like -> DESI-like (secondary axis = primary * yfactor)
JUMP    = dict(ymin=5e-10, ylow=5e-9, yhigh=1e-4, ymax=5, gap=0.5)   # broken y axis, same for both panels

#one entry per panel; the file name AND the printed R_min come from the same numbers (no way to get out of sync)
panels = [dict(mmin=0.5, rmin=49.8), dict(mmin=1.0, rmin=29.3)]

#columns of the input files
COL = dict(sigma=0, abh_full=1, bns_full=2, pbh_full=3, pbh_full_err=4,
           abh_hod=5, abh_obs=6, bns_hod=7, bns_obs=8, pbh_frac=9, pbh_frac_err=10)

#------------------------------------------------------------------ helpers
def add_converted_yaxis(ax, factor, cut=None, labels=True):
	"""Secondary y axis on the right showing  y*factor,  with exactly the same scale as `ax`
	(log, or broken log from add_y_jump).  Call it AFTER add_y_jump / set_yscale on `ax`.
	cut=(ylow,yhigh) is the removed range in the units of `ax`: no ticks are put inside it.

	Why not secondary_yaxis(functions=...)+set_yscale('log'): matplotlib's own secondary scale only
	knows the linear conversion, and set_yscale('log') then spreads the converted range as a plain
	log axis, which does not match the jump scale of the parent -> ticks end up in the wrong place."""
	tr  = ax.yaxis.get_transform()          # parent data -> scaled coordinate (includes the jump)
	itr = tr.inverted()
	to_scaled   = lambda s: tr.transform_non_affine(np.asarray(s, dtype=float)/factor)
	from_scaled = lambda v: itr.transform_non_affine(np.asarray(v, dtype=float))*factor
	sec = ax.secondary_yaxis('right', functions=(lambda y: y*factor, lambda y: y/factor))
	sec.set_yscale('function', functions=(to_scaled, from_scaled))   # secondary value s sits where parent value s/factor sits
	#ticks: full decades of the converted quantity, only in the visible part
	lo, hi  = np.array(ax.get_ylim())*factor
	ticks   = 10.0**np.arange(np.ceil(np.log10(lo)), np.floor(np.log10(hi))+1)
	if cut is not None:
		ticks = ticks[(ticks/factor <= cut[0]) | (ticks/factor >= cut[1])]
	sec.yaxis.set_major_locator(mticker.FixedLocator(ticks))
	sec.yaxis.set_major_formatter(mticker.LogFormatterSciNotation())
	sec.yaxis.set_minor_locator(mticker.NullLocator())
	sec.tick_params(axis='y', which='both', labelright=labels)
	return sec


def plot_events(ax, ev, sigma):
	"""All curves of one panel. ev = events array (columns as in COL), already multiplied by `years`."""
	fade = dict(marker='.', fade_dex=0.3, alpha0=0.3)
	s12  = slice(1, 3)    #ABH only has the two middle sigma values
	#ABH
	plot_upper_limit(ax, sigma[s12], ev[COL['abh_full']][s12], 'blue', r'ABH full voids',     **fade)
	plot_upper_limit(ax, sigma[s12], ev[COL['abh_hod']][s12],  'blue', r'ABH HOD void frac',  ls='--', **fade)
	plot_upper_limit(ax, sigma[s12], ev[COL['abh_obs']][s12],  'blue', r'ABH Obs void frac',  ls=':',  **fade)
	#BNS
	plot_upper_limit(ax, sigma, ev[COL['bns_full']], 'green', r'BNS full voids',    **fade)
	plot_upper_limit(ax, sigma, ev[COL['bns_hod']],  'green', r'BNS HOD void frac', ls='--', **fade)
	plot_upper_limit(ax, sigma, ev[COL['bns_obs']],  'green', r'BNS Obs void frac', ls=':',  **fade)
	#PBH
	ax.plot(sigma, ev[COL['pbh_frac']], 'r--', label=r'PBH void frac')
	ax.errorbar(sigma, ev[COL['pbh_frac']], yerr=ev[COL['pbh_frac_err']], fmt='r.')
	ax.plot(sigma, ev[COL['pbh_full']], 'r-', label=r'PBH full voids')
	ax.errorbar(sigma, ev[COL['pbh_full']], yerr=ev[COL['pbh_full_err']], fmt='r.')


#------------------------------------------------------------------ figure
fig, axs = plt.subplots(1, len(panels), figsize=(9, 5), sharey=True, layout='constrained')
axs[0].set_ylabel('BOSS-like\n\\# event / 10yr', fontsize=14)    #newline instead of '\\' (a LaTeX '\\' in a label is fragile)

for ax, p in zip(axs, panels):
	ev = np.loadtxt(f"sigma_Nabh_Nbns_Npbh_Mmin{p['mmin']:.1f}_Rmin{p['rmin']}_{config}.txt").T
	ev[1:] *= years
	plot_events(ax, ev, ev[COL['sigma']])

	ax.set_xlabel(r'$\sigma$', fontsize=16)
	#text in axes coordinates (independent of the broken y scale)
	ax.text(0.97, 0.16, r"$M_{\rm min} = %g\,M_\odot$" % p['mmin'], transform=ax.transAxes, fontsize=13, ha='right', va='top')
	ax.text(0.97, 0.08, r"$R_{\rm min} = %g$ Mpc" % p['rmin'],      transform=ax.transAxes, fontsize=13, ha='right', va='top')

#y jump: the axes are shared, so scale, limits and tick locators are common; the cut-out mask and the // marks
#are per-axes artists, hence one call per panel
for ax in axs:
	add_y_jump(ax, **JUMP)

#secondary (DESI-like) axes: ticks on both panels, labels only on the last one
cut = (JUMP['ylow'], JUMP['yhigh'])
for i, ax in enumerate(axs):
	sec = add_converted_yaxis(ax, yfactor, cut=cut, labels=(i == len(axs)-1))
	if i == len(axs)-1:
		sec.set_ylabel('DESI-like\n\\# event / 10yr', fontsize=14, rotation=270, labelpad=28)

#single legend for both panels, title carries the configuration
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='outside upper center', ncols=3, frameon=False,
           prop={'size': 12}, title="ET configuration 2L", title_fontsize=18)

plt.savefig(f"Number_of_events_{config}.pdf", bbox_inches='tight')
plt.show()
