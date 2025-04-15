import numpy as np
import cmocean

yCrossSection = 3000
xCrossSection = 8000
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
saltRange = np.linspace(31, 35, 128)
tempRange = np.linspace(-2, 3.0, 128)
uRange = np.linspace(-.2, .2, 127)
vRange = np.linspace(-.2, .2, 127)
wRange = np.linspace(-.005, .005, 127)
meltRange = np.linspace(0,.5,128)