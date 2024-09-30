from MITgcmutils import mds
import imageio
import numpy as np
import matplotlib.pyplot as plt

saveGif = True

XC = mds.rdmds('XC')
YC = mds.rdmds('YC')
ZC = mds.rdmds('RC')
figDir = 'figs/'
ZC = ZC.squeeze()
XC = XC.squeeze()
YC = YC.squeeze()
print('XC:',np.shape(XC))
print('YC:',np.shape(YC))
print('ZC:',np.shape(ZC))

startStep = 864
endStep = 7776
stepSize = 864

for di in [2,4,6,8,10,12,14]:
    zSlice = di
    depth = ZC[zSlice]
    print(depth)

    for i in np.arange(stepSize, endStep+1, stepSize):
        print(i)
        T = mds.rdmds('T', i)
        #print('T:',np.shape(T))
        T = T[zSlice,:,:].squeeze()
        #print('T:',np.shape(T))
        #print("T: ", np.min(T) , " " , np.max(T))

        W = mds.rdmds('W', i)
        #print('W:',np.shape(W))
        W = W[zSlice,:,:].squeeze()
        #print('W:',np.shape(W))
        #print("W: ", np.min(W) , " " , np.max(W))

        U = mds.rdmds('U', i)
        U = U[zSlice,:,:].squeeze()
        #print("U: ", np.min(U) , " " , np.max(U))

        S = mds.rdmds('S', i)
        S = S[zSlice,:,:].squeeze()
        #print("S: ", np.min(S) , " " , np.max(S))

        plt.figure
        plt.contourf(XC, YC, T, np.linspace(-1.2, 3, 64), cmap='hot_r')
        plt.colorbar()
        plt.title('Temperature @ %f step %05i' % (depth,i))
        plt.savefig("%sT%05i.png" % (figDir, i))
        # plt.show()
        plt.close()

        plt.figure
        plt.contourf(XC, YC, U, np.linspace(-0.05, 0.05, 64), cmap='bwr')
        plt.colorbar()
        plt.title('U velocity @ %f step %05i' % (depth,i))
        plt.savefig("%sU%05i.png" % (figDir, i))
        # plt.show()
        plt.close()

        plt.figure
        plt.contourf(XC, YC, W, np.linspace(-1e-3, 1e-3, 64), cmap='bwr')
        plt.colorbar()
        plt.title('W velocity @ % step %05i' % (depth,i))
        plt.savefig("%sW%05i.png" % (figDir, i))
        # plt.show()
        plt.close()

        # plt.figure(3)
        # plt.quiver(XC, ZC, U, W)
        # plt.colorbar()
        # plt.title('Flow Quiver Plot')
        # plt.gca().set_aspect('equal')
        # plt.show()

        plt.figure(4)
        plt.contourf(XC, YC, S, np.linspace(33, 35, 64), cmap='bwr')
        plt.colorbar()
        plt.title('Salt @ % step %05i' % (depth,i))
        plt.savefig("%sS%05i.png" % (figDir, i))
        # plt.show()
        plt.close()

    if(saveGif):
        toPlot = ['U','W','T','S']
        for q in toPlot:
            filenames = ['']
            images = []
            for i in np.arange(stepSize, endStep+1, stepSize):
                if i == stepSize:
                    filenames = ['figs/%s%05i.png' % (q,i)]
                else:
                    filenames.append('figs/%s%05i.png' % (q,i))

            print(filenames)

            for filename in filenames:
                images.append(imageio.imread(filename))
            imageio.mimsave('figs/%03i%s.gif' % (np.abs(depth), q), images)
        