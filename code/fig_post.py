#! /usr/bin/env python

import numpy as np
import sys
import os
import scipy.interpolate as si
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.ioff()

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.utils as mu
sys.path.append("../pyratbay/modules/MCcubed/MCcubed/plots")
import colormaps as cm

# Set palette color:
palette = cm.viridis_r
palette.set_under(color='w')
palette.set_bad(color='w')

rjup = 7.1492e4
burnin  = 5000
nbins   = 35
nlevels = 20

# Fixed-Na MCMC run:
post = np.load("../run/MCMC_WASP17b_hires.npz")

Z = post["Z"]
Z = Z[burnin:,:]
Z[:,1] /= rjup
nsamples, npars1 = np.shape(Z)
posterior1 = Z[Z[:,0]>0,:]
# Gather 2D histograms:
hist1 = []
xran1, yran1, lmax1 = [], [], []
for   j in np.arange(1, npars1): # Rows
  for i in np.arange(npars1-1):  # Columns
    if j > i:
      h,x,y = np.histogram2d(posterior1[:,i], posterior1[:,j],
                             bins=nbins, normed=False)
      hist1.append(h.T)
      xran1.append(x)
      yran1.append(y)
      lmax1.append(np.amax(h)+1)

# Free-Na MCMC run:
post = np.load("../run/MCMC_WASP17b_hires_abundance.npz")
Z = post["Z"]
Z = Z[burnin:,:]
Z[:,1] /= rjup
nsamples, npars2 = np.shape(Z)
posterior2 = Z[Z[:,0]>0,:]
# Gather 2D histograms:
hist2 = []
xran2, yran2, lmax2 = [], [], []
for   j in np.arange(1, npars2): # Rows
  for i in np.arange(npars2-1):  # Columns
    if j > i:
      h,x,y = np.histogram2d(posterior2[:,i], posterior2[:,j],
                             bins=nbins, normed=False)
      hist2.append(h.T)
      xran2.append(x)
      yran2.append(y)
      lmax2.append(np.amax(h)+1)


# Plotting setup:
pdf  = [None]*npars2
xpdf = [None]*npars2
percentile = 0.683
hkw = {'histtype':'step', 'lw':1.5}
parname = [r"$T\ ({\rm K})$", r"$R_0\ (R_{\rm Jup})$", r"[Na/H]"]
fs = 10


fig = plt.figure(21, (8.5, 4))
plt.clf()
plt.subplots_adjust(0.095, 0.17, 0.98, 0.95, wspace=0.1, hspace=0.1)

# Fixed-Na MCMC run:
xran = [(450, 2050),(1.74, 1.88)]
xticks = [(500, 1000, 1500, 2000), (1.75, 1.8, 1.85)]
# Pairwise posterior:
ax=plt.subplot(256)
k = 0
a = plt.contourf(hist1[k], cmap=palette, vmin=1, origin='lower',
            levels=[0]+list(np.linspace(1,lmax1[k], nlevels)),
            extent=(xran1[k][0], xran1[k][-1], yran1[k][0], yran1[k][-1]))
plt.xticks(size=fs-1, rotation=90)
plt.yticks(size=fs-1)
plt.ylabel(parname[1], size=fs, multialignment='center')
plt.xlabel(parname[0], size=fs)
ax.set_xlim(xran[0])
ax.set_ylim(xran[1])
ax.set_xticks(xticks[0])
for c in a.collections:
  c.set_edgecolor("face")
# Marginal posteriors:
for i in np.arange(npars1):
  ax = plt.subplot(2, 5, 1+i*6)
  a = plt.xticks(size=fs-1, rotation=90)
  a = plt.yticks(size=fs-1)
  ax.set_yticks([])
  if i == 0:
    ax.set_xticklabels("")
  else:
    plt.xlabel(parname[i], size=fs)
  vals, bins, h = plt.hist(posterior1[:,i], bins=25, normed=False, **hkw)
  if percentile is not None:
    PDF, Xpdf, HPDmin = mu.credregion(posterior1[:,i], percentile,
                                      pdf[i], xpdf[i])
    vals = np.r_[0, vals, 0]
    bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
    f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
    ax.fill_between(Xpdf, 0, f(Xpdf), where=PDF>=HPDmin,
                 facecolor='0.7', edgecolor='none', interpolate=False)
  ax.set_xlim(xran[i])
  ax.set_xticks(xticks[i])

# Free-Na MCMC run:
xran = [(350,3000), (1.68, 1.95), (-4.05, 3.05)]
# Pairwise posterior:
idx = 10, 16, 17
for k in np.arange(npars2):
  ax = plt.subplot(3, 6, idx[k])
  a = plt.contourf(hist2[k], cmap=palette, vmin=1, origin='lower',
            levels=[0]+list(np.linspace(1,lmax2[k], nlevels)),
            extent=(xran2[k][0], xran2[k][-1], yran2[k][0], yran2[k][-1]))
  plt.yticks(size=fs-1)
  plt.xticks(size=fs-1, rotation=90)
  if k <2:
    plt.ylabel(parname[k+1], size=fs, multialignment='center')
  else:
    ax.set_yticklabels("")
  if k > 0:
    plt.xlabel(parname[k-1], size=fs)
  else:
    ax.set_xticklabels("")
  for c in a.collections:
    c.set_edgecolor("face")
# Marginal posteriors:
for i in np.arange(npars2):
  ax = plt.subplot(3, 6, 4+i*7)
  a = plt.xticks(size=fs-1, rotation=90)
  a = plt.yticks(size=fs-1)
  ax.set_yticks([])
  if i == npars2-1:
    plt.xlabel(parname[i], size=fs)
  else:
    ax.set_xticklabels("")
  vals, bins, h = plt.hist(posterior2[:,i], bins=25, normed=False, **hkw)
  if percentile is not None:
    PDF, Xpdf, HPDmin = mu.credregion(posterior2[:,i], percentile,
                                      pdf[i], xpdf[i])
    vals = np.r_[0, vals, 0]
    bins = np.r_[bins[0] - (bins[1]-bins[0]), bins]
    f = si.interp1d(bins+0.5*(bins[1]-bins[0]), vals, kind='nearest')
    ax.fill_between(Xpdf, 0, f(Xpdf), where=PDF>=HPDmin,
                    facecolor='0.7', edgecolor='none', interpolate=False)
for j in np.arange(npars2):
  for i in np.arange(j+1):
    ax = plt.subplot(3, 6, 4+i+6*j)
    ax.set_xlim(xran[i])
    if j>i:
      ax.set_ylim(xran[j])
# Colorbar:
bounds = np.linspace(0, 1.0, nlevels)
norm = mpl.colors.BoundaryNorm(bounds, palette.N)
ax2 = fig.add_axes([0.9, 0.525, 0.015, 0.4])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=palette, norm=norm,
      spacing='proportional', boundaries=bounds, format='%.1f')
cb.set_label("Posterior density", fontsize=fs)
cb.set_ticks(np.linspace(0, 1, 5))
for c in ax2.collections:
  c.set_edgecolor("face")

plt.savefig("../plots/WASP-17b_hires_posteriors.pdf")



