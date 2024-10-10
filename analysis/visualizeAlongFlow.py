from MITgcmutils import mds
import imageio
import numpy as np
import matplotlib.pyplot as plt
import os

dataDir = 'shelf500/'
saveGif = True
figDir = '%sfigs/' % dataDir
dataDir = '%sresults/' % dataDir
print(figDir)
print(dataDir)

x = mds.rdmds('%sXC' % dataDir)
y = mds.rdmds('%sYC' % dataDir)
z = mds.rdmds('%sRC' % dataDir)

maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir(dataDir):
    # print(file)
    if "T.0" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])
        if int(words[1]) < startStep and int(words[1]) > 0:
            sizeStep = int(words[1])
            startStep = int(words[1])

if(maxStep/sizeStep > 50):   #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale
print(startStep,sizeStep,maxStep)

os.system('rm figs/sideDump*.png')
os.system('rm figs/autoSideDump*.gif')

crossSection = 1250

ySlice = np.argmin(np.abs(y[:,0] - crossSection))
print('cross section is y =', y[ySlice,0])
width = y[ySlice,0]

x = x[ySlice,:].squeeze()
z = z.squeeze()
for i in np.arange(startStep, maxStep+1, sizeStep):
    print(i)
    T = mds.rdmds('%sT' % dataDir, i)
    #print('T:',np.shape(T))
    T = T[:,ySlice,:].squeeze()
    #print('T:',np.shape(T))
    print("T: ", np.min(T) , " " , np.max(T))

    W = mds.rdmds('%sW' % dataDir, i)
    #print('W:',np.shape(W))
    W = W[:,ySlice,:].squeeze()
    #print('W:',np.shape(W))
    print("W: ", np.min(W) , " " , np.max(W))

    U = mds.rdmds('%sU' % dataDir, i)
    U = U[:,ySlice,:].squeeze()
    print("U: ", np.min(U) , " " , np.max(U))

    S = mds.rdmds('%sS' % dataDir, i)
    S = S[:,ySlice,:].squeeze()
    print("S: ", np.min(S) , " " , np.max(S))

    V = mds.rdmds('%sV' % dataDir, i)
    V = V[:,ySlice,:].squeeze()
    print("V: ", np.min(V) , " " , np.max(V))

    plt.figure
    plt.contourf(x, z, T, np.linspace(-1, 1, 128), cmap='hot_r')
    plt.colorbar()
    plt.title('Temperature @ %f step %05i' % (width,i))
    plt.savefig("%ssideTdump%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

    plt.figure
    plt.contourf(x, z, U, np.linspace(-0.01, 0.01, 128), cmap='bwr')
    plt.colorbar()
    plt.title('U velocity @ %f step %05i' % (width,i))
    plt.savefig("%ssideUdump%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

    plt.figure
    plt.contourf(x, z, W, np.linspace(-0.0001, 0.0001, 128), cmap='bwr')
    plt.colorbar()
    plt.title('W velocity @ % step %05i' % (width,i))
    plt.savefig("%ssideWdump%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

    plt.figure
    plt.contourf(x, z, V, np.linspace(-6e-5, 6e-5, 128), cmap='bwr')
    plt.colorbar()
    plt.title('W velocity @ % step %05i' % (width,i))
    plt.savefig("%ssideVdump%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

    # plt.figure(3)
    # plt.quiver(XC, ZC, U, W)
    # plt.colorbar()
    # plt.title('Flow Quiver Plot')
    # plt.gca().set_aspect('equal')
    # plt.show()

    plt.figure(4)
    plt.contourf(x, z, S, np.linspace(33, 35, 128), cmap='viridis')
    plt.colorbar()
    plt.title('Salt @ % step %05i' % (width,i))
    plt.savefig("%ssideSdump%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

if(saveGif):
    toPlot = ['U','W','T','S','V']
    for q in toPlot:
        os.system('magick -delay 50 %sside%s*.png -colors 256 -depth 256 %sautoSide%03i%s.gif' %(figDir,q,figDir,np.abs(width),q))
        print('magick %ssidedump%s*.png %ssidedump%03i%s.gif' %(figDir,q,figDir,np.abs(width),q))
        # filenames = ['']
        # images = []
        # for i in np.arange(startStep, endStep+1, stepSize):
        #     if i == startStep:
        #         filenames = ['%sside%s%05i.png' % (figDir,q,i)]
        #     else:
        #         filenames.append('%sside%s%05i.png' % (figDir,q,i))

        # print(filenames)

        # for filename in filenames:
        #     images.append(imageio.v2.imread(filename))
        # imageio.mimsave('%sside%03i%s.gif' % (figDir,np.abs(width), q), images, palettesize=128, duration=500)
    