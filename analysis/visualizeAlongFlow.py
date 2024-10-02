from MITgcmutils import mds
import imageio
import numpy as np
import matplotlib.pyplot as plt
import os

dataDir = 'B95_1000_864/'
saveGif = True
figDir = '%sfigs/' % dataDir
dataDir = '%sresults/' % dataDir
print(figDir)
print(dataDir)

XC = mds.rdmds('%sXC' % dataDir)
YC = mds.rdmds('%sYC' % dataDir)
ZC = mds.rdmds('%sRC' % dataDir)

ZC = ZC.squeeze()
XC = XC.squeeze()
YC = YC.squeeze()
print('XC:',np.shape(XC))
print('YC:',np.shape(YC))
print('ZC:',np.shape(ZC))

startStep = 00
endStep = 864
stepSize = 864
ySlice = 3
XC = XC[ySlice,:].squeeze()
width = YC[ySlice,0]
print(width)

for i in np.arange(startStep, endStep+1, stepSize):
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

    plt.figure
    plt.contourf(XC, ZC, T, np.linspace(-1.2, 3, 128), cmap='hot_r')
    plt.colorbar()
    plt.title('Temperature @ %f step %05i' % (width,i))
    plt.savefig("%ssideT%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

    plt.figure
    plt.contourf(XC, ZC, U, np.linspace(-0.007, 0.007, 128), cmap='bwr')
    plt.colorbar()
    plt.title('U velocity @ %f step %05i' % (width,i))
    plt.savefig("%ssideU%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

    plt.figure
    plt.contourf(XC, ZC, W, np.linspace(-6e-5, 6e-5, 128), cmap='bwr')
    plt.colorbar()
    plt.title('W velocity @ % step %05i' % (width,i))
    plt.savefig("%ssideW%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

    # plt.figure(3)
    # plt.quiver(XC, ZC, U, W)
    # plt.colorbar()
    # plt.title('Flow Quiver Plot')
    # plt.gca().set_aspect('equal')
    # plt.show()

    plt.figure(4)
    plt.contourf(XC, ZC, S, np.linspace(33, 35, 128), cmap='bwr')
    plt.colorbar()
    plt.title('Salt @ % step %05i' % (width,i))
    plt.savefig("%ssideS%05i.png" % (figDir, i))
    # plt.show()
    plt.close()

if(saveGif):
    toPlot = ['U','W','T','S']
    for q in toPlot:
        os.system('magick -delay 50 %sside%s*.png -colors 256 -depth 256 %sside%03i%s.gif' %(figDir,q,figDir,np.abs(width),q))
        print('magick %sside%s*.png %sside%03i%s.gif' %(figDir,q,figDir,np.abs(width),q))
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
    