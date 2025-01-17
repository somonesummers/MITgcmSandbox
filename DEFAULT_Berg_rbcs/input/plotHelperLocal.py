import numpy as np
import cmocean

yCrossSection = 1500
xCrossSection = 15000
zDepth = -50
plotDPI = 175
usePcolor = False
showZeros = True

# Color Maps
saltCmap = "cmo.haline"
tempCmap = "cmo.thermal"
uCmap = "cmo.balance"
vCmap = "cmo.balance"
wCmap = "cmo.curl"
meltCmap = "cmo.rain"

# Color Ranges
saltRange = np.linspace(32, 35, 128)
tempRange = np.linspace(-9, 0, 128)
uRange = np.linspace(-.5, .5, 127)
vRange = np.linspace(-.5, .5, 127)
wRange = np.linspace(-.05, .05, 127)
meltRange = np.linspace(0,.5,128)