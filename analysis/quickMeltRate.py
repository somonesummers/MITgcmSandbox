# %%
from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import cmocean

ddir = ''

os.system('rm -f figs/fwf*.png')
os.system('rm -f figs/autoFwf*.gif')

maxStep = 0
sizeStep = 1e10
startStep = 1e10


if(os.path.isfile('input/data.iceberg')):
    isShelf = False
    loadFile = "BRGFlx"
    name = ['BRGfwFlx','BRGhtFlx','BRGmltRt'] 
    legendName = ['Fresh Water Flux [kg/m^s/s]','Heat Flux [W/m^2]','MeltRate [m/yr]'] 
else:
    isShelf = True
    loadFile = "SHIflux"
    name = ['SHIfwFlx','SHIhtFlx','SHIForcT','SHIForcS'] 
    legendName = ['Fresh Water Flux [m/yr]','Heat Flux [W/m^2]','SHIForcT','SHIForcS'] 
     

for file in os.listdir('%sresults' % ddir):
    # print(file)
    if loadFile in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])
        if int(words[1]) < startStep and int(words[1]) > 0:
            sizeStep = int(words[1])
            startStep = int(words[1])

if(maxStep/sizeStep > 50):  #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print(startStep,sizeStep,maxStep)

# %%
#---------+----+----+----------------+-----------------
# <-Name->|Levs|grid|<--  Units   -->|<- Tile (max=80c)
#---------+----+----+----------------+-----------------
# SHIfwFlx|  1 |SM  |kg/m^2/s        |Ice shelf fresh water flux (positive upward)
# SHIhtFlx|  1 |SM  |W/m^2           |Ice shelf heat flux  (positive upward)
# SHIUDrag| 30 |UU  |m/s^2           |U momentum tendency from ice shelf drag
# SHIVDrag| 30 |VV  |m/s^2           |V momentum tendency from ice shelf drag
# SHIForcT|  1 |SM  |W/m^2           |Ice shelf forcing for theta, >0 increases theta
# SHIForcS|  1 |SM  |g/m^2/s         |Ice shelf forcing for salt, >0 increases salt

# %%
y = mds.rdmds("%sresults/YC" % ddir)
x = mds.rdmds("%sresults/XC" % ddir)
z = mds.rdmds("%sresults/RC" % ddir)


topo = np.fromfile('input/topog.slope', dtype='>f8')
# ice = np.fromfile('input/icetopo.exp1', dtype='>f8')
topo = topo.reshape(np.shape(x))
# ice = ice.reshape(np.shape(x))

# iceEdge = np.interp(z[0,0,0],ice[0,:],x[0,:])


for i in np.arange(startStep, maxStep + 1, sizeStep):
    data = mds.rdmds("%sresults/%s" % (ddir, loadFile ), i)
    if(isShelf):
        dataTemp = np.squeeze(data[0, :, :]*1e-3*3.154e7) #kg/(m^2*s) -> m/yr
        levels = np.linspace(-70, 0, 128)
        cMap = "cmo.rain_r"
    else:
        dataTemp = np.max(np.squeeze(data[2,:,:,:]),0)
        levels = np.linspace(0, 1, 128)
        cMap = cm = "cmo.rain"
    for k in [0]:
        lvl = levels
        cm = cMap
        cp = plt.contourf(
            np.squeeze(x),
            np.squeeze(y),
            dataTemp,
            lvl,
            extend="both",
            cmap=cm
        )
    # plt.plot(iceEdge*np.ones(np.shape(x)),np.squeeze(y[:,0]),color='gray')
    cbar = plt.colorbar(cp)
    cbar.set_label('[m/yr]')
    plt.xlabel('Along Fjord [m] %.2f %.2f nan: %i' %(np.nanmin(dataTemp),np.nanmax(dataTemp),np.max(np.isnan(dataTemp))))
    plt.ylabel('Across Fjord [m]')
    plt.title("Max Meltrate at %i" % (i))
    j = i/startStep
        
    str = "%sfigs/fwf%05i.png" % (ddir,j)
        
    plt.savefig(str, format='png')
    
    #plt.show()
    plt.close()

os.system('magick -delay 5 %sfigs/fwf*.png -colors 256 -depth 256 %sfigs/autoFwf.gif' %(ddir, ddir))

# %%


# %%



