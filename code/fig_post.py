#! /usr/bin/env python

import numpy as np
import sys
import os
import scipy.interpolate as si
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.ioff()

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.plots as mp
import MCcubed.utils as mu
sys.path.append("../pyratbay/modules/MCcubed/MCcubed/plots")
import colormaps as cm

# Set palette color:
palette = cm.viridis_r
palette.set_under(color='w')
palette.set_bad(color='w')

rjup = 7.1492e4
burnin = 5000
#post = np.load("../retrieval_v03_3HW/MCMC_WASP17b_hires.npz")
post = np.load("../run/MCMC_WASP17b_hires.npz")
Z = post["Z"]
Z = Z[burnin:,:]
Z[:,1] /= rjup

parname = [r"$T\ ({\rm K})$", r"$R_0\ (R_{\rm Jup})$"]
nsamples, npars = np.shape(Z)
thinning = 1
nbins = 35
nlevels = 20
posterior = Z[:-1]

# Gather 2D histograms:
hist = []
xran, yran, lmax = [], [], []
for   j in np.arange(1, npars): # Rows
  for i in np.arange(npars-1):  # Columns
    if j > i:
      h,x,y = np.histogram2d(posterior[0::thinning,i],
                  posterior[0::thinning,j], bins=nbins, normed=False)
      hist.append(h.T)
      xran.append(x)
      yran.append(y)
      lmax.append(np.amax(h)+1)


# Histogram keywords depending whether one wants the HPD or not:
pdf  = [None]*npars
xpdf = [None]*npars
hkw = {}
percentile = 0.683
hkw = {'histtype':'step', 'lw':2}
ranges = [None]*npars

fs = 11
fig = plt.figure(20, (8.5, 3))
plt.clf()
plt.subplots_adjust(0.095, 0.25, 0.98, 0.95, wspace=0.4)
ax=plt.subplot(131)

h = 1 # Subplot index
k = 0 # Histogram index
plt.yticks(size=fs)
plt.ylabel(parname[1], size=fs, multialignment='center')
plt.xticks(size=fs, rotation=90)
plt.xlabel(parname[0], size=fs)
# The plot:
a = plt.contourf(hist[0], cmap=palette, vmin=1, origin='lower',
            levels=[0]+list(np.linspace(1,lmax[k], nlevels)),
            extent=(xran[k][0], xran[k][-1], yran[k][0], yran[k][-1]))
for c in a.collections:
  c.set_edgecolor("face")
# The colorbar:
bounds = np.linspace(0, 1.0, nlevels)
norm = mpl.colors.BoundaryNorm(bounds, palette.N)
ax2 = fig.add_axes([0.34, 0.25, 0.015, 0.7])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=palette, norm=norm,
      spacing='proportional', boundaries=bounds, format='%.1f')
cb.set_label("Posterior density", fontsize=fs)
cb.set_ticks(np.linspace(0, 1, 5))
for c in ax2.collections:
  c.set_edgecolor("face")

for i in np.arange(npars):
  ax = plt.subplot(1, 3, i+2)
  a = plt.xticks(size=fs-2.0, rotation=90)
  a = plt.yticks(size=fs-2.0)
  ax.set_yticklabels([])
  plt.ylabel(r"$N\ \rm samples$", fontsize=fs)
  plt.xlabel(parname[i], size=fs)
  vals, bins, h = plt.hist(posterior[0::thinning, i], bins=25,
                           range=ranges[i], normed=False, **hkw)
  # Plot HPD region:
  if percentile is not None:
    PDF, Xpdf, HPDmin = mu.credregion(posterior[:,i], percentile,
                                      pdf[i], xpdf[i])
    vals = np.r_[0, vals, 0]
    bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
    # interpolate xpdf into the histogram:
    f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
    # Plot the HPD region as shaded areas:
    if ranges[i] is not None:
      Xran = np.argwhere((Xpdf>ranges[i][0]) & (Xpdf<ranges[i][1]))
      Xpdf = Xpdf[np.amin(Xran):np.amax(Xran)]
      PDF  = PDF [np.amin(Xran):np.amax(Xran)]
    ax.fill_between(Xpdf, 0, f(Xpdf), where=PDF>=HPDmin,
                 facecolor='0.7', edgecolor='none', interpolate=False)
plt.savefig("../plots/posteriors_v02.ps")


