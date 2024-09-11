from MITgcmutils import mds
import numpy as np
import matplotlib.pyplot as plt
XC = mds.rdmds('XC')
YC = mds.rdmds('YC')
ZC = mds.rdmds('RC')

ZC = ZC.squeeze()
XC = XC.squeeze()

endStep = 864000

T = mds.rdmds('T', endStep)
T = T[:,0,:].squeeze()

W = mds.rdmds('W', endStep)
W = W[:,0,:].squeeze()
print("W: ", np.min(W) , " " , np.max(W))

U = mds.rdmds('U', endStep)
U = U[:,0,:].squeeze()
print("U: ", np.min(U) , " " , np.max(U))

S = mds.rdmds('S', endStep)
S = S[:,0,:].squeeze()
print("S: ", np.min(S) , " " , np.max(S))

plt.figure(1)
plt.contourf(XC, ZC, T, np.linspace(-1.8, .4, 12), cmap='hot_r')
plt.colorbar()
plt.title('Temperature')
plt.savefig("temp.png")
plt.show()

plt.figure(2)
plt.contourf(XC, ZC, U, np.linspace(-0.12, 0.12, 12), cmap='bwr')
plt.colorbar()
plt.title('U velocity')
plt.savefig("U.png")
plt.show()

plt.figure(3)
plt.contourf(XC, ZC, W, np.linspace(-0.04, 0.04, 12), cmap='bwr')
plt.colorbar()
plt.title('W velocity')
plt.savefig("W.png")
plt.show()

# plt.figure(3)
# plt.quiver(XC, ZC, U, W)
# plt.colorbar()
# plt.title('Flow Quiver Plot')
# plt.gca().set_aspect('equal')
# plt.show()

# plt.figure(4)
# plt.contourf(XC, ZC, S, np.linspace(33, 35, 64), cmap='bwr')
# plt.colorbar()
# plt.title('Salt?')
# plt.show()