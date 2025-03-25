import numpy as np
import cmocean
from matplotlib.colors import LinearSegmentedColormap

def red_yellow_white_cyan_blue():
    N = 256
    cols = ['#000055', '#0000f5', '#008cff',
            # '#7affff', '#ffffff', '#ffff83',
            '#9ee6e6', '#ececec', '#eaea99',
            '#ff9d0c', '#ff0500', '#5f0000']
    return LinearSegmentedColormap.from_list('custom', cols, N)


yCrossSection = 1500
xCrossSection = 15000
zDepth = -50
plotDPI = 175
usePcolor = False
showZeros = True

# Color Maps
saltCmap = "cmo.haline"
tempCmap = "cmo.thermal"
uCmap = red_yellow_white_cyan_blue()
vCmap = "cmo.balance"
wCmap = "cmo.curl"
meltCmap = "cmo.rain"


# Color Ranges
saltRange = np.linspace(32, 35, 128)
tempRange = np.linspace(-9, 0, 128)
uRange = np.linspace(-.15, .15, 127)
vRange = np.linspace(-.5, .5, 127)
wRange = np.linspace(-.05, .05, 127)
meltRange = np.linspace(0,.5,128)