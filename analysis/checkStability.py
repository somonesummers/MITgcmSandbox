import fileinput
import numpy as np

u_char = .5 #[m/s]
deltaT = 0
f0 = 0
viscAz = 0
dx = 500 #[m]
dz = 6.66
nz = 75

for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            deltaT = float(line[8:-2])
        elif "viscAz=" in line:
            viscAz = float(line[8:-2])
        elif "f0" in line:
            f0 = float(line[4:-2])
        elif "f0" in line:
            f0 = float(line[4:-2])

openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
iceCoverage = 1 - np.nanmin(openFrac[openFrac > 0])

print('deltaT', deltaT, ', viscAz', viscAz, ', f0', f0, ', dx', dx,', iceCoverage',iceCoverage)







u_char = .5 #[m/s]
S_adv = 2 * (u_char * deltaT)/dx  # < 0.5 (Courant–Friedrichs–Lewy)
S_in = f0 * deltaT # < 0.5 (adams bashforth II)
# S_lh = 8 * params01['viscAh'] * deltaT /(run_config['horiz_res_m']**2) # < 0.6 
S_lv = 4 * viscAz * deltaT /(dz**2) # < 0.6 

S_adv_brg = 2 * (u_char * deltaT)/(dx*(1 - iceCoverage))  # < 0.5 (Courant–Friedrichs–Lewy)
S_in = f0 * deltaT # < 0.5 (adams bashforth II)
S_lv_brg = 4 * viscAz * deltaT /((dz*(1 - iceCoverage))**2) # < 0.6
print('====== Stability Check =====')
print("S_adv: <0.5 , S_in: <0.5 , S_lv: <0.6")
print("S_adv: %.04f, S_in: %.04f, S_lv: %.06f" % (S_adv, S_in, S_lv))
print("S_adv: %.04f, S_in: %.04f, S_lv: %.06f" % (S_adv_brg, S_in, S_lv_brg))