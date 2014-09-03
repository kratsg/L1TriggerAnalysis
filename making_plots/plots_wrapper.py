import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm

#for getting FWHM
from scipy.interpolate import UnivariateSpline, interp1d

#for getting error function fitting in diff/int curves
from scipy.optimize import curve_fit
from scipy.special import erf, erfinv

from matplotlib.ticker import FuncFormatter

class PlotHelpers(object):
  def __init__(self, **kwargs):
    self.dataSetStr = kwargs.get('dataSetStr', None)
    self.towerThrStr = kwargs.get('towerThrStr', None)
    self.seedCutStr = kwargs.get('seedCutStr', None)
    self.noiseCutStr = kwargs.get('noiseCutStr', None)
    self.figsize = kwargs.get('figsize', (16, 12) )
    self.labelsize = kwargs.get('labelsize', 30)
    self.titlesize = kwargs.get('titlesize', 36)
    self.ticksize = kwargs.get('ticksize', 26)
    self.linewidth = kwargs.get('linewidth', 4)
    self.light_grey = np.array([float(200)/float(255)]*3)
    self.cmap = kwargs.get('cmap', pl.cm.autumn)
    self.textprops = dict(boxstyle='round', facecolor=self.light_grey, alpha=0.5, linewidth=0.0)
    self.regions = {'all': [-4.9,  4.9],\
                        1: [-1.6,  0.0],\
                        2: [ 0.0,  1.6],\
                        3: [-4.9, -1.6],\
                        4: [ 1.6,  4.9],\
                     '3a': [-2.5, -1.6],\
                     '3b': [-3.2, -2.5],\
                     '3c': [-4.9, -3.2],\
                     '4a': [ 1.6,  2.5],\
                     '4b': [ 2.5,  3.2],\
                     '4c': [ 3.2,  4.9]}

    self.labels = {'rho.gFEX': r'gFEX $\rho$ [GeV/area]',\
                   'rho.offline': r'Offline $\rho$ [GeV/area]',\
                   'vxp.number': r'Number of primary vertices (vxp_nTracks $\geq$ 2)',\
                   'gTower.Et': r'$E_T^\mathrm{gTower}$ [GeV]',\
                   'gTower.multiplicity': 'gTower multiplicity / event'}

  def region_legend(self, region):
    return r'towers: $\eta \in [{:0.1f}, {:0.1f})$'.format(*self.regions[region])

  def binomial_errors(hist_ratio, hist_one, hist_two):
    errors = []
    for w, num, den in zip(hist_ratio, hist_one, hist_two):
      # root.cern.ch/root/html/src/TH1.cxx.html#l5.yxD
      # formula cited (for histograms [num, den] with no errors) is:
      #     w = num/den
      #     if w = 1:
      #             sigma = 0
      #     else:
      #             sigma = abs( (1 - 2*w + w**2) / den**2 )
      if w == 1.0:
        errors.append(0.0)
      else:
        errors.append( (np.abs( (1.-2.*w + w**2.)/den**2.))**0.5/2. )
    return errors

  def _profile(self, xbins, xvals, yvals):
    # this finds the difference between successive bins and then halves it, and
    #     then adds it back to the bins to get xpoints
    # gets the average yvalue for each of the xbins
    xpoints = xbins[:-1] + np.array([v-xbins[i-1] for i,v in enumerate(xbins)][1:])/2.
    # this digitizes our samples by figuring out which xbin each xval belongs to
    digitized = np.digitize(xvals, xbins)
    ymean = [np.mean(yvals[np.where(digitized == i)]) for i in np.arange(1,len(xbins)+1)][1:]
    #filter out missing values
    nonnan = np.where(~np.isnan(ymean))
    xpoints = np.array(xpoints)[nonnan]
    ymean = np.array(ymean)[nonnan]
    return xpoints, ymean

  def profile_y(self, xbins, xvals, yvals):
    return self._profile(xbins, xvals, yvals)

  def profile_x(self, ybins, yvals, xvals):
    return self._profile(ybins, yvals, xvals)
   
  def FWHM(bins, vals):
    spline = UnivariateSpline(bins[:-1]+np.diff(bins)/2., vals-np.max(vals)/2., s=0)
    roots = spline.roots() # find the roots
    r1, r2 = roots[0], roots[-1]
    return np.abs(r1-r2)

  #this is a wrapper around file strings to ensure the directory exists
  def write_file(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
      os.makedirs(d)
    return f

  def btwn(self, val, a, b):
    return (a <= val)&(val < b)

  def region_cut(self, data, index, region):
    if region == 1:
      return self.btwn(data[index], -1.6, 0.0)
    elif region == 2:
      return self.btwn(data[index], 0.0, 1.6)
    elif region == 3:
      return self.btwn(data[index], -4.9, -1.6)
    elif region == 4:
      return self.btwn(data[index], 1.6, 4.9)
    elif region == '3a':
      return self.btwn(data[index], -2.5, -1.6)
    elif region == '3b':
      return self.btwn(data[index], -3.2, -2.5)
    elif region == '3c':
      return self.btwn(data[index], -4.9, -3.2)
    elif region == '4a':
      return self.btwn(data[index], 1.6, 2.5)
    elif region == '4b':
      return self.btwn(data[index], 2.5, 3.2)
    elif region == '4c':
      return self.btwn(data[index], 3.2, 4.9)

  def add_legend(self, fig, axes):
    legend = ax.legend(fancybox=True, framealpha=0.75, fontsize=self.labelsize)
    legend.get_frame().set_facecolor(light_grey)
    legend.get_frame().set_linewidth(0.0)

  def add_description(self, fig, ax, va='top', ha='left', strings=[]):
    if va=='top':
      ypos = 0.95
    elif va=='bottom':
      ypos = 0.05
    if ha=='left':
      xpos = 0.05
    elif ha=='right':
      xpos = 0.95
    ax.text(xpos, ypos, "\n".join(strings), transform=ax.transAxes, fontsize=self.labelsize, verticalalignment=va, horizontalalignment=ha, bbox=self.textprops)

  def add_atlas(self, fig, ax, level=0):
    if level == 0:
      textstr = 'Internal'
    elif level == 1:
      textstr = 'Preliminary'
    ax.text(0.05, 0.95, 'ATLAS', fontsize=42, style='italic', fontweight='bold', verticalalignment='top', horizontalalignment='left', transform=fig.gca().transAxes)
    ax.text(0.27, 0.95, textstr, verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)
    ax.text(0.05, 0.90, 'Simulation', verticalalignment='top', horizontalalignment='left', fontsize=40, transform=fig.gca().transAxes)
    return fig, ax


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
    ax.grid(True, which='both', linewidth=3, linestyle='--', alpha=0.5)

  def corr2d(self, x, y, bins_x, bins_y, label_x, label_y, profile_x=False, profile_y=False):
    corr = np.corrcoef(x, y)[0,1]
    fig, ax = pl.subplots(figsize=self.figsize)
    hist, _, _, _ = ax.hist2d(x, y, norm=LogNorm(), bins=(bins_x, bins_y) , alpha=0.75, cmap = self.cmap)

    cbar = ax.colorbar(ticks=np.logspace(0, np.log10(np.max(hist)), 10), format=self.label_formatter)
    cbar.set_label('number density', fontsize=self.labelsize)
    cbar.ax.tick_params(labelsize=self.ticksize)

    ax.xlabel(label_x, fontsize=self.labelsize)
    ax.ylabel(label_y, fontsize=self.labelsize)

    self.add_grid(fig, ax)

    corrStr = '$\mathrm{Corr} = {:0.4f}f$'.format(corr)
    self.add_description(fig, ax, va='bottom', ha='right', strings=[self.dataSetStr, self.towerThrStr, corrStr])

    self.add_atlas(fig, ax, level=1)

    ax.tick_params(axis='both', which='major', labelsize=self.ticksize)
    return fig, ax
    pl.savefig( write_file('ZH/%s_%s_correlation.png' % (filename_id, 'gFEX_rho_1')), bbox_inches='tight')
    pl.savefig( write_file('ZH/%s_%s_correlation.eps' % (filename_id, 'gFEX_rho_1')), bbox_inches='tight')
    pl.savefig( write_file('ZH/%s_%s_correlation.pdf' % (filename_id, 'gFEX_rho_1')), bbox_inches='tight')
    pl.close()

  def turnon(self, den, num):
    pass

