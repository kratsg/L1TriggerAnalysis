import os

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm

# for getting FWHM
from scipy.interpolate import UnivariateSpline

# for getting error function fitting in diff/int curves
from scipy.optimize import curve_fit
from scipy.special import erf, erfinv

from matplotlib.ticker import FuncFormatter


class PlotHelpers(object):
  def __init__(self, **kwargs):
    self.dataSetStr = kwargs.get('dataSetStr', None)
    self.towerThrStr = kwargs.get('towerThrStr', None)
    self.seedCutStr = kwargs.get('seedCutStr', None)
    self.noiseCutStr = kwargs.get('noiseCutStr', None)
    self.figsize = kwargs.get('figsize', (16, 12))
    self.labelsize = kwargs.get('labelsize', 30)
    self.titlesize = kwargs.get('titlesize', 36)
    self.ticksize = kwargs.get('ticksize', 26)
    self.linewidth = kwargs.get('linewidth', 4)
    self.light_grey = np.array([float(200)/float(255)]*3)
    self.cmap = kwargs.get('cmap', pl.cm.hot_r)
    self.textprops = dict(boxstyle='round', facecolor=self.light_grey, alpha=0.5, linewidth=0.0)
    self.regions = {'all': [-4.9,  4.9],
                    1: [-1.6,  0.0],
                    2: [+0.0,  1.6],
                    3: [-4.9, -1.6],
                    4: [+1.6,  4.9],
                    '3a': [-2.5, -1.6],
                    '3b': [-3.2, -2.5],
                    '3c': [-4.9, -3.2],
                    '4a': [+1.6,  2.5],
                    '4b': [+2.5,  3.2],
                    '4c': [+3.2,  4.9]}
    self.colors = ['b', 'r', 'c', 'm', 'y', 'k']
    self.markers = ['o', 'v', '^', 's', 'D', '8']

    self.labels = {'rho.gFEX': r'gFEX $\rho$ [GeV/area]',
                   'rho.offline': r'Offline $\rho$ [GeV/area]',
                   'vxp.number': r'Number of primary vertices (vxp_nTracks $\geq$ 2)',
                   'gTower.Et': r'$E_T^\mathrm{gTower}$ [GeV]',
                   'gTower.multiplicity': 'gTower multiplicity / event'}

  def btwn(self, val, a, b):
    if b is None:
      return (a <= val)
    elif a is None:
      return (val < b)
    else:
      return (a <= val) & (val < b)

  def region_cut(self, data, region):
    return self.btwn(data, *self.regions[region])

  def offline_region_cut(self, data, region):
    return self.btwn(data, -1.0 + self.regions[region][0], 1.0+self.regions[region][1])

  def region_legend(self, region):
    return r'towers: $\eta \in [{:0.1f}, {:0.1f})$'.format(*self.regions[region])

  def binomial_errors(self, hist_ratio, hist_one, hist_two):
    errors = []
    for w, num, den in zip(hist_ratio, hist_one, hist_two):
      # root.cern.ch/root/html/src/TH1.cxx.html#l5.yxD
      # formula cited (for histograms [num, den] with no errors) is:
      #     e1, e2 are errors of the bins, either by sqrt(sum(weights**2 for
      #     bin)) or poisson sqrt(bin)
      #     w = num/den
      #     if w = 1:
      #             sigma = 0
      #     else:
      #             sigma = abs( ((1 - 2*w)* e1**2 + w**2 * e2**2) / den**2 )
      if w == 1.0:
        errors.append(0.0)
      else:
        errors.append((np.abs(((1.-2.*w)*num + (w**2.)*den)/den**2.))**0.5)
    return errors

  def _profile(self, xbins, xvals, yvals):
    # this finds the difference between successive bins and then halves it, and
    #     then adds it back to the bins to get xpoints
    # gets the average yvalue for each of the xbins
    shift = np.array([v-xbins[i-1] for i, v in enumerate(xbins)])/2.
    # assume that the shift is constant
    shift[0] = shift[-1]
    # determine the xpoints by shifting
    xpoints = xbins + shift
    # this digitizes our samples by figuring out which xbin each xval belongs to
    digitized = np.digitize(xvals, xbins)  # bins by bins[i-1] <= x < bins[i], left inclusive, right exclusive
    num_yvals = np.array([len(yvals[np.where(digitized == i)]) for i in np.arange(1, len(xbins)+1)])
    # filter out missing values by num_yvals > 0
    xpoints = xpoints[np.where(num_yvals > 0)]
    ymean = np.array([np.mean(yvals[np.where(digitized == i)]) for i in np.arange(1, len(xbins)+1)[np.where(num_yvals > 0)]])
    yerr  = np.array([np.std(yvals[np.where(digitized == i)]) for i in np.arange(1, len(xbins)+1)[np.where(num_yvals > 0)]])
    return xpoints, ymean, yerr

  def profile_y(self, xbins, xvals, yvals):
    return self._profile(xbins, xvals, yvals)

  def profile_x(self, ybins, yvals, xvals):
    return self._profile(ybins, yvals, xvals)

  def FWHM(self, bins, vals):
    spline = UnivariateSpline(bins[:-1]+np.diff(bins)/2., vals-np.max(vals)/2., s=0)
    roots = spline.roots()  # find the roots
    r1, r2 = roots[0], roots[-1]
    return np.abs(r1-r2)

  # this is a wrapper around file strings to ensure the directory exists
  def write_file(self, f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
      os.makedirs(d)
    return f

  def add_legend(self, fig, ax, **kwargs):
    legend = ax.legend(fancybox=True, framealpha=0.75, fontsize=self.labelsize, **kwargs)
    if isinstance(legend, None.__class__): return False
    legend.get_frame().set_facecolor(self.light_grey)
    legend.get_frame().set_linewidth(0.0)

  def add_description(self, fig, ax, align='tl', strings=[]):
    if align[0] == 't':
      ypos = 0.95
    elif align[0] == 'b':
      ypos = 0.05
    elif align[0] == 'c':
      ypos = 0.45
    if align[1] == 'l':
      xpos = 0.05
    elif align[1] == 'r':
      xpos = 0.95

    va = {'t': 'top', 'b': 'bottom', 'c': 'center'}[align[0]]
    ha = {'l': 'left', 'r': 'right'}[align[1]]

    ax.text(xpos, ypos, "\n".join(filter(None, strings)), transform=ax.transAxes, fontsize=self.labelsize, verticalalignment=va, horizontalalignment=ha, bbox=self.textprops)

  def add_atlas(self, fig, ax, level=0):
    if level == 0:
      textstr = 'Internal'
    elif level == 1:
      textstr = 'Preliminary'
    ax.text(0.05, 0.95, 'ATLAS', fontsize=42, style='italic', fontweight='bold', verticalalignment='top', horizontalalignment='left', transform=fig.gca().transAxes)
    ax.text(0.27, 0.95, textstr, verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)
    ax.text(0.05, 0.90, 'Simulation', verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)

  def latex_float(self, f, pos=0):
      float_str = "{0:.2g}".format(f)
      if "e" in float_str:
          base, exponent = float_str.split("e")
          return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
      else:
          return r"${}$".format(float_str)

  @property
  def label_formatter(self):
    return FuncFormatter(self.latex_float)

  def add_grid(self, fig, ax):
    ax.grid(True, which='major', linewidth=3, linestyle='--', alpha=0.5)

  def add_labels(self, fig, ax, xlabel=None, ylabel=None, title=None):
    if xlabel is not None:
      ax.set_xlabel(xlabel, fontsize=self.labelsize)
    if ylabel is not None:
      ax.set_ylabel(ylabel, fontsize=self.labelsize)
    if title is not None:
      ax.set_title(title, fontsize=self.titlesize)
    ax.tick_params(axis='both', which='both', labelsize=self.ticksize)

  def format_cbar(self, cbar, label='number density'):
    cbar.set_label(label, fontsize=self.labelsize, labelpad=40)
    cbar.ax.tick_params(labelsize=self.ticksize)

  def to_file(self, fig, ax, filename, transparent=True):
    fig.savefig(self.write_file(filename), bbox_inches='tight', transparent=transparent)

  def corr2d(self, x, y, bins_x, bins_y, label_x, label_y, xlim=None, ylim=None, profile_x=False, profile_y=False, title=None, strings=[], align='bl', ticks=None, weights=None):
    if weights is None:
      weights = np.ones(x.size)
    corr = np.corrcoef(x, y)[0, 1]
    fig, ax = pl.subplots(figsize=self.figsize)
    counts, edges_x, edges_y, im = ax.hist2d(x, y, bins=(bins_x, bins_y), norm=LogNorm(), alpha=0.75, cmap = self.cmap, weights=weights)
    ticks = ticks or np.logspace(0, np.log10(np.max(counts)), 10)
    cbar = fig.colorbar(im, ticks=ticks, format=self.label_formatter)
    self.format_cbar(cbar)

    corrStr = '$\mathrm{{Corr}} = {:0.4f}$'.format(corr)
    self.add_description(fig, ax, align=align, strings=strings+[corrStr])

    self.add_labels(fig, ax, xlabel=label_x, ylabel=label_y, title=title)

    if profile_y:
      points_x, mean_y, err_y = self.profile_y(edges_x, x, y)
      ax.scatter(points_x, mean_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)
    if profile_x:
      points_y, mean_x, err_x = self.profile_x(edges_y, y, x)
      ax.scatter(mean_x, points_y, s=80, facecolor='w', edgecolor='k', marker='o', linewidth=2)

    self.add_grid(fig, ax)

    if xlim is not None:
      ax.set_xlim(xlim)
    if ylim is not None:
      ax.set_ylim(ylim)

    return fig, ax

  def add_turnon(self, fig, ax, data=None, den_cut=None, num_cut=None, label=None, bins=np.arange(0.0, 1000.0, 10.0), kind='differential', p0=(0., 0., 0., 0.)):

    # *_noPileup means the data with the pileup subtracted
    hist_efficiency_den, _ = np.histogram(data[np.where(den_cut)], bins=bins)
    hist_efficiency_num, _ = np.histogram(data[np.where(den_cut & num_cut)], bins=bins)

    nonzero_bins = np.where(hist_efficiency_den != 0)
    # compute integral and differential curves
    if kind == 'differential':
      denominator = hist_efficiency_den[nonzero_bins]
      numerator   = hist_efficiency_num[nonzero_bins]
    else:
      denominator = np.cumsum(hist_efficiency_den[nonzero_bins][::-1])[::-1]
      numerator   = np.cumsum(hist_efficiency_num[nonzero_bins][::-1])[::-1]

    hist_eff_curve = np.true_divide(numerator, denominator)

    # halfway between bins is where we plot
    width_efficiency = np.array([x - bins[i-1] for i, x in enumerate(bins)][1:])
    xpoints_efficiency = bins[:-1] + width_efficiency/2.

    # binomial errors s^2 = n * p * q
    errors_eff = self.binomial_errors(hist_eff_curve, numerator, denominator)

    try:
      w = self.fit_func(xpoints_efficiency[nonzero_bins], hist_eff_curve, p0=p0)
    except TypeError:  # fix later
      w = (-1,-1,-1,-1)
    ax.errorbar(xpoints_efficiency[nonzero_bins], hist_eff_curve, yerr=errors_eff, ecolor='black', label='{}\n{}'.format(label, self.make_turnon_label(w)), linewidth=self.linewidth)
    return xpoints_efficiency, hist_eff_curve, errors_eff, nonzero_bins, w

  def erf_func(self, x, a, b, c, d):
    # note that b == sigma here, see wiki for more info
    return a*erf((x-c)/b) + d

  def fit_func(self, x, y, p0=(0., 0., 0., 0.)):
    try:
      # define erfx used for error fitting
      popt, pcov = curve_fit(self.erf_func, x, y, p0=p0)
    except RuntimeError:
      return -1.
    return popt  # return (a,b,c,d); width = b

  def peak_point(self, w):
    return w[1]*erfinv((0.95-w[3])/w[0])+w[2]

  def make_turnon_label(self, w):
    if np.all(w == -1) or ~np.isfinite(self.peak_point(w)):
      return ''
    else:
      return '$w = {0:0.1f},\ x_{{0.95}} = {1:0.1f}$'.format(w[1], self.peak_point(w))


# start defining a meta-class to register all instances of classes made
# http://jakevdp.github.io/blog/2012/12/01/a-primer-on-python-metaclasses/#Example-2:-Registering-Subclasses
class Plottable(type):
  # we use __init__ rather than __new__ here because we want to
  # modify attributes of the class *after* they have been created
  def __init__(cls, name, bases, dct):
    if not hasattr(cls, 'registry'):
      # this is the base class. Create an empty registry
      cls.registry = {}
    else:
      if cls.registerMe:
        # this is a derived class. Add cls to the registry
        interface_id = name.lower()
        cls.registry[interface_id] = cls

    super(Plottable, cls).__init__(name, bases, dct)


class L1TriggerAnalysis(object):
  __metaclass__ = Plottable
  registerMe = False

  @classmethod
  def draw(cls, event):
    raise NotImplementedError("{}.draw() should be implemented".format(cls.__name__))

  @classmethod
  def region_draw(cls, event):
    raise NotImplementedError("{}.region_draw() should be implemented".format(cls.__name__))
