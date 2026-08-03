"""
Custom diverging colormap matching the paper's own Fig. 3 heatmaps (GBEK,
PRB 88, 235106 (2013)), sampled directly from the published colorbar
(page_09.png) rather than approximated from a standard matplotlib map --
red (high) -> orange -> yellow -> green -> blue -> white (~0) -> purple
(low), not a plain red-white-blue diverging scale like RdBu_r.

GBEK_CMAP is defined over the paper's own value range [-0.2, 0.5] (its
Fig. 3 colorbar limits); GBEK_VMIN/GBEK_VMAX are exposed alongside it so
callers reproduce the exact same visual mapping the paper uses by default,
but can pass their own vmin/vmax to imshow(cmap=GBEK_CMAP, ...) if their
data's natural range differs.
"""
from matplotlib.colors import LinearSegmentedColormap

GBEK_VMIN = -0.2
GBEK_VMAX = 0.5

_CONTROL_POINTS = [
    (0.00000, (0.6549, 0.4745, 0.7647)),
    (0.03129, (0.7059, 0.5333, 0.7961)),
    (0.06257, (0.7529, 0.5882, 0.8196)),
    (0.09371, (0.8000, 0.6471, 0.8471)),
    (0.12500, (0.8353, 0.7020, 0.8706)),
    (0.15629, (0.8784, 0.7608, 0.8980)),
    (0.18757, (0.9176, 0.8196, 0.9216)),
    (0.21871, (0.9529, 0.8824, 0.9490)),
    (0.25000, (0.9804, 0.9490, 0.9765)),
    (0.28129, (0.9961, 0.9961, 1.0000)),
    (0.31257, (0.8471, 0.8157, 0.9255)),
    (0.34371, (0.6667, 0.6235, 0.8431)),
    (0.37500, (0.4431, 0.4275, 0.7451)),
    (0.40629, (0.3216, 0.3294, 0.6824)),
    (0.43757, (0.2588, 0.3020, 0.6627)),
    (0.46871, (0.2314, 0.2902, 0.6510)),
    (0.50000, (0.2078, 0.2627, 0.6314)),
    (0.53129, (0.1608, 0.2471, 0.6039)),
    (0.56257, (0.1294, 0.3137, 0.4549)),
    (0.59371, (0.1294, 0.4235, 0.2706)),
    (0.62500, (0.1569, 0.5412, 0.2353)),
    (0.65629, (0.1882, 0.6431, 0.2510)),
    (0.68743, (0.2549, 0.7137, 0.2353)),
    (0.71871, (0.3725, 0.7608, 0.1961)),
    (0.75000, (0.5294, 0.8078, 0.1490)),
    (0.78129, (0.6941, 0.8588, 0.0941)),
    (0.81257, (0.8706, 0.9098, 0.0588)),
    (0.84371, (0.9725, 0.9137, 0.0275)),
    (0.87500, (0.9804, 0.7059, 0.0392)),
    (0.90629, (0.9647, 0.5059, 0.0745)),
    (0.93743, (0.9451, 0.3176, 0.1059)),
    (0.96871, (0.9373, 0.1922, 0.1294)),
    (1.00000, (0.9294, 0.1255, 0.1373)),
]

GBEK_CMAP = LinearSegmentedColormap.from_list("gbek_paper", _CONTROL_POINTS, N=256)
