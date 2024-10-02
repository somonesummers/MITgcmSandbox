# This is a python script that generates the input data
# variable x resolution

import numpy as np
import sys
import matplotlib.pyplot as plt


def write_bin(fname, data):
    print(fname, np.shape(data))
    data.astype(">f8").tofile(fname)


# Dimensions of grid
H = 200.0
Lx = 25.0e3
Ly = 5.0e3

nx = 50
ny = 10
nz = 40

dx = Lx / nx
dy = Ly / ny
dz = H / nz
print('dx:',dx," dy:",dy," dz:",dz)
# Some constants
gravity = 9.81
sbeta = 8.0e-4
talpha = 0.4e-4
rho0 = 999.8
T0 = 1
S0 = 35

x = np.zeros([ny, nx])
x[:, 0] = dx / 2

for i in np.arange(1, nx):
    x[:,i] = x[:, i - 1] + dx

z = -dz / 2 + np.arange(0, -H, -dz)

# Temperature profile
tcd = 50
Tmin = -1.6
Tmax = 0.2
Tc = (Tmax + Tmin) / 2
Trange = Tmax - Tmin
T2 = np.zeros([nz,ny])
for j in np.arange(0,ny):
    T2[:,j] = Tc - Trange / 2 * np.tanh(np.pi * (z + tcd) / tcd)
Tconst = np.zeros([nz,ny]) + 0.4
T = T2

Sc = 34.5
Srange = -1
S2 = np.zeros([nz,ny])
for j in np.arange(0,ny):
    S2[:,j] = Sc + Srange / 2 * np.tanh(np.pi * (z + tcd) / tcd)
Sconst = np.zeros([nz,ny]) + 35

S = S2

Rref = rho0 * (1 - talpha * (T - T0) + sbeta * (S - S0))

t = np.zeros([nz, ny, nx])
s = np.zeros([nz, ny, nx])

for j in np.arange(0,ny):
    for k in np.arange(0, nz):
        t[k, j, :] = t[k, j, :] + T[k,j]
        s[k, j, :] = s[k, j, :] + S[k,j]

ubound = np.zeros([nz, ny])

write_bin("T.init", t)
write_bin("S.init", s)
write_bin("S.bound", S)
write_bin("T.bound", T)
write_bin("U.bound", ubound)


plt.plot(S - 34, z, 'b', label="Sref - 34")
plt.plot(T, z, 'r', label="Tref")
plt.legend()
plt.savefig("initialTS")
plt.show()

# Topographie

d = np.zeros([ny, nx]) - H
# d[:, 0] = 0
# d[:, -1] = 0
write_bin("topog.slope", d)

# Ice shelf
gldepth = -H + 10
icebase = -10
sgdd = 0
L = [0, 100, 15000]
D = [-H + sgdd, gldepth, icebase]

m = np.zeros([ny, nx]) - 1
iceshelf = np.zeros([ny, nx])
icefront = np.zeros([ny, nx])
m = np.zeros([ny, nx])

for j in np.arange(0, ny):
    # print(j)
    for i in np.arange(1, np.shape(D)[0]):
        xi = np.where((x[j,:] >= L[i - 1]) & (x[j,:] < L[i]))[0]
        m[j, xi] = (D[i] - D[i - 1]) / (L[i] - L[i - 1])
        iceshelf[j, xi] = np.linspace(D[i - 1], D[i], np.size(xi))
        # print(np.linspace(D[i - 1], D[i], np.size(xi)))
        icefront[j, xi] = np.linspace(D[i - 1], D[i], np.size(xi))

iceshelf[:, 0] = D[0]
icefront[:, 0] = D[0]

iceshelf[x > L[-1]] = 0
icefront[x > L[-1]] = 0
# iceshelf[m > dz / dx] = np.nanmin(iceshelf[m < dz / dx])
icefront[m < dz / dx] = 0
icelength = np.zeros([ny, nx])
icelength[icefront < 0] = 1 / dy

a0 = -0.0573
c0 = 0.0832
b = -7.53e-4

if sgdd != 0:
    sgdu = np.zeros([ny, nz])
    sgdS = np.zeros([ny, nz])
    sgdT = np.zeros([ny, nz])
    sgd = 0.001
    sgdi = np.where((z < gldepth) & (z >= gldepth - sgdd))
    sgdu[:, sgdi] = sgd / dy / sgdd
    sgdS[:, sgdi] = 30
    sgdT[:, sgdi] = a0 * sgdS[:, sgdi] + c0 + b * -z[sgdi]
    write_bin("T.sgd", sgdT)
    write_bin("S.sgd", sgdS)
    write_bin("U.sgd", sgdu)


plt.plot(np.transpose(x), np.transpose(iceshelf), 'r', label="shelfice")
# plot(sgdu * 1000, z)
# plt.plot(x, icefront[:, 0], label="icefront")
plt.plot(np.transpose(x), np.transpose(d), 'b', label="bathy")
plt.legend()
plt.savefig("geo")
plt.show()

write_bin("icetopo.exp1", iceshelf)
write_bin("icefrontdepth", icefront)
write_bin("icefrontlength", icelength)

# Phi 0

pano = np.zeros([ny, nx])
for j in np.arange(0,ny):
    for i in np.arange(0, nx):
        ki = np.where(z >= iceshelf[j,i])[0]

        if not ki.size > 0:
            pextra = 0
            panoex = 0
            ptop = 0
        else:
            k = np.nanmax(ki)
            ptop = np.sum(Rref[0:k,j] * gravity * dz)  # Ice pressure
            ptopano = ptop - rho0 * gravity * dz * (k + 1)  # Ice pressure anomaly
            pextra = abs(z[k] - iceshelf[j,i]) * gravity * rho0
            panoex = pextra - abs(z[k] - iceshelf[j,i]) * gravity * Rref[k,j]

        pano[j, i] = panoex + ptopano


write_bin("phi0.exp1", pano)
