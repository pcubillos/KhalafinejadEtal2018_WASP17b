#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ioff()

sys.path.append("../pyratbay/modules/MCcubed")
import MCcubed.utils as mu

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants as pc


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# High resolution RT model:
pyrat = pb.pbay.run("spectrum_bestfit.cfg")

vel = pyrat.phy.rvsystem
wn  = pyrat.spec.wn / ((1-vel/pc.c) / np.sqrt(1-(vel/pc.c)**2))
bandwl = 1/(pyrat.obs.bandwn[1:]*pc.A)
uncert = pyrat.obs.uncert[1:]
data   = pyrat.obs.data[1:]

# N-sigma boundaries around Na lines:
N = 3.0
wn0 = np.array(pyrat.alkali.model[0].wn)
wn1 = wn0 / ((1-vel/pc.c) / np.sqrt(1-(vel/pc.c)**2))
wl1 = 1/(wn1*pc.A)
dwl = 1/(wn0*pc.A) / pyrat.phy.resolution
wlran1 = wl1[0] - N*dwl[0], wl1[0] + N*dwl[0]
wlran2 = wl1[1] - N*dwl[1], wl1[1] + N*dwl[1]


# Best-fit value spectrum:
pyrat.atm.temp[:] = 1550.0
pyrat.phy.rplanet = 129200 * pc.km
pyrat = pb.pyrat.run(pyrat, [pyrat.atm.temp, pyrat.atm.q])
# Double normalization of high-resolution data:
pyrat.spec.hires = gaussf((1-pyrat.spec.spectrum) / (1-pyrat.spec.spectrum[0]),
                          pyrat.phy.sigma)
pyrat.obs.bandflux = si.interp1d(wn,pyrat.spec.hires)(pyrat.obs.bandwn)
flux = np.copy(pyrat.obs.bandflux[1:])

# Spectra for R0 5% larger/smaller:
pyrat.phy.rplanet = 1.05 * 129200 * pc.km
pyrat = pb.pyrat.run(pyrat, [pyrat.atm.temp, pyrat.atm.q])
pyrat.spec.hires = gaussf((1-pyrat.spec.spectrum) / (1-pyrat.spec.spectrum[0]),
                          pyrat.phy.sigma)
pyrat.obs.bandflux = si.interp1d(wn,pyrat.spec.hires)(pyrat.obs.bandwn)
hiflux = np.copy(pyrat.obs.bandflux[1:])

pyrat.phy.rplanet = 0.95 * 129200 * pc.km
pyrat = pb.pyrat.run(pyrat, [pyrat.atm.temp, pyrat.atm.q])
pyrat.spec.hires = gaussf((1-pyrat.spec.spectrum) / (1-pyrat.spec.spectrum[0]),
                          pyrat.phy.sigma)
pyrat.obs.bandflux = si.interp1d(wn,pyrat.spec.hires)(pyrat.obs.bandwn)
loflux = np.copy(pyrat.obs.bandflux[1:])


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# The figure:
fs = 12
lw = 1.5
plt.figure(3, (8.5, 3.1))
plt.clf()
plt.subplots_adjust(0.1, 0.15, 0.95, 0.97, wspace=0.05)
ax=plt.subplot(121)
plt.plot(bandwl, flux,   color="b",         lw=lw)
plt.plot(bandwl, hiflux, color="orangered",    lw=lw, zorder=-5)
plt.plot(bandwl, loflux, color="forestgreen", lw=lw, zorder=-5)
plt.errorbar(bandwl, data, uncert, fmt="ko", ms=3, ecolor="0.3")
plt.axvspan(wlran2[0], wlran2[1], color="0.85", zorder=-6)
xticks = [5887, 5888, 5889, 5890, 5891]
ax.set_xticks(xticks)
ax.set_xticklabels(np.asarray(xticks, str))
plt.xlim(5887.25, 5890.75)
plt.ylim(0.97, 1.02)
plt.ylabel(r"$R$", fontsize=fs)
plt.ylabel(r"$\tilde{\mathfrak{R}}$", fontsize=fs)
plt.xlabel(r"$\rm Wavelength\ \ (A)$", fontsize=fs)
# Right panel:
ax=plt.subplot(122)
plt.plot(bandwl, flux,   color="b",         lw=lw)
plt.plot(bandwl, hiflux, color="orangered", lw=lw, zorder=-5)
plt.plot(bandwl, loflux, color="forestgreen", lw=lw, zorder=-5)
plt.errorbar(bandwl, data, uncert, fmt="ko", ms=3, ecolor="0.3")
plt.axvspan(wlran1[0], wlran1[1], color="0.85", zorder=-6)
ax.set_yticklabels([])
xticks = [5893, 5894, 5895, 5896, 5897]
ax.set_xticks(xticks)
ax.set_xticklabels(np.asarray(xticks, str))
plt.xlim(5893.25, 5896.75)
plt.ylim(0.97, 1.02)
plt.xlabel(r"$\rm Wavelength\ \ (A)$", fontsize=fs)
plt.savefig("../plots/WASP-17b_hires_fix-Na_spectrum.pdf")
