__author__ = 'Chris Jones'

"""
Reusable no burnout equation page.

Units used:
gram
centimeter
second
"""

import numpy as np
from scipy.integrate import quad
import sys
import matplotlib.pyplot as plt
sys.path.insert(0, 'C:/Users/Chris Jones/Documents/PythonProjects/OFHC_Therm_Cond')
import heat_xfer as hx

I0 = 200.  # Initial current in amps
If = 1.  # Final current in amps
u0 = 293.  # Temp at start of dump for hottest portion of conductor in K
uf = 493.  # Temp desired after slow dump for hottest portion of conductor inK
deltaT = uf - u0  # Allowable temperature shift during dump
tauFast = 15.0  # Fast dump time constant
tauSlow = 150.0  # Slow dump time constant
r273 = 1.71e-6  # Conductor resistivity at 293K in Ohm * cm
rrr = 150  # The ratio of resistivity from 273 K to 4 K
c500 = 0.408  # Specific Heat of copper at 500K, most conservative in deltaT
alpha = 3.86e-3  # Copper's temperature coefficient of resistance per Kelvin
rho = 8.96  # Density of copper in g/cm^3
Lcoil = 10.  # Coil inductance in Henries

# Selectable analyses
r_c_constant = 0
slow_dump = 0
fast_dump = 0
cw = 1
temp_var = 1
cw_and_temp = 0

###
# Calculate values used in further analysis
###
# Resistances needed to achieve given time constant tau
rSlow = Lcoil / tauSlow
rFast = Lcoil / tauFast
# print('rSlow = {:.2e} Ohm\nrFast = {:.2e} Ohm'.format(rSlow, rFast))

# Time needed to go from I0 to If
tSlow = tauSlow * np.log(I0 / If)
tFast = tauFast * np.log(I0 / If)
# print('tSlow = {:.2f} s\ntFast = {:.2f} s'.format(tSlow, tFast))


# Function for finding area given conservative constant resistivity and specific heat in a slow dump
if r_c_constant:
    c493 = hx.specific_heat(493)
    A0 = I0 * np.sqrt((alpha * r273 * tauSlow)
                      / (2 * rho * c493 * np.log((alpha * uf + 1)/(alpha * u0 + 1))))
    print('\nUsing conservative, constant specific heat no burnout area is: {:.4f} cm^2'.format(A0))
    print('For a conductor diameter of {:.4f} cm'.format(np.sqrt(4 * A0 / np.pi)))

# Using temp dependent specific heat and resistivity for slow dump
if slow_dump:
    tempNumerator = lambda t: np.exp(-2 * t / tauSlow)
    numerator = quad(tempNumerator, 0, tSlow)

    tempDenominator = lambda u: hx.specific_heat(u) / hx.resistivity_BG(u, rrr=rrr, rho273=r273)
    denominator = quad(tempDenominator, u0, uf)

    A = I0 * np.sqrt(numerator[0] / (rho * denominator[0]))

    print('\nUsing temp dependent specific heat and resistivity, minimum area is {:.4f} cm^2'.format(A))
    print('For a minimum conductor diameter of {:.4f} cm'.format(np.sqrt(4 * A / np.pi)))

# Using temp dependent specific heat and resistivity for fast dump
if fast_dump:
    tempNumerator = lambda t: np.exp(-2 * t / tauFast)
    numerator = quad(tempNumerator, 0, tFast)

    tempDenominator = lambda u: hx.specific_heat(u) / hx.resistivity_BG(u, rrr=rrr, rho273=r273)
    denominator = quad(tempDenominator, u0, uf)

    A = I0 * np.sqrt(numerator[0] / (rho * denominator[0]))

    print('\nUsing temp dependent specific heat and resistivity, minimum area is {:.4f} cm^2'.format(A))
    print('For a minimum conductor diameter of {:.4f} cm'.format(np.sqrt(4 * A / np.pi)))

# Analyzing Cold Work's effect on resistivity and burnout area
if cw:
    delta_r = np.linspace(0., 0.6e-7, num=25)
    # delta_r = np.array([0., 0.3e-7, 0.45e-7])
    A = []

    for i in range(len(delta_r)):
        tempNumerator = lambda t: np.exp(-2 * t / tauSlow)
        numerator = quad(tempNumerator, 0, tSlow)

        tempDenominator = lambda u: hx.specific_heat(u) / hx.resistivity_BG(u, rrr=rrr, rho273=r273, delta_r=delta_r[i])
        denominator = quad(tempDenominator, u0, uf)

        A.append(I0 * np.sqrt(numerator[0] / (rho * denominator[0])))

    A = np.array(A)
    print(A)
    b = np.sqrt(4 * A / np.pi)
    print(b)
    print(A[-1]-A[0])
    fig = plt.figure()
    ax = fig.gca()
    ax.annotate(r'25% CW',
                    arrowprops=dict(facecolor='black', shrink=0.02, width=2),
                    xytext=(2e-8, 0.1030), textcoords='data', size=20,
                    xy=(0.3e-7, 0.1035), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    ax.annotate(r'50% CW',
                    arrowprops=dict(facecolor='black', shrink=0.02, width=2),
                    xytext=(4e-8, 0.1032), textcoords='data', size=20,
                    xy=(0.45e-7, 0.10382), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    plt.plot(delta_r, A, 'k')
    plt.grid(True)
    plt.xlabel(r'Increase in Resistivity $(\Omega cm)$')
    plt.ylabel(r'No Burnout Area ($cm^2$)')
    plt.title(r'Cold Working Effect on No Burnout Minimum Area')
    plt.show()

if temp_var:
    delta_u = np.linspace(280, 305, num=25)
    A = np.zeros(len(delta_u))

    tempNumerator = lambda t: np.exp(-2 * t / tauSlow)
    numerator = quad(tempNumerator, 0, tSlow)

    for i in range(len(delta_u)):
        tempDenominator = lambda u: hx.specific_heat(u) / hx.resistivity_BG(u, rrr=rrr, rho273=r273)
        denominator = quad(tempDenominator, delta_u[i], uf)

        A[i] = I0 * np.sqrt(numerator[0] / (rho * denominator[0]))

    print(A)

    plt.plot(delta_u, A, 'k')
    plt.grid(True)
    plt.xlabel(r'Temperature (K)')
    plt.ylabel(r'No Burnout Area ($cm^2$)')
    plt.title(r'Temperature Variation Effect on No Burnout Minimum Area')
    plt.show()

if cw_and_temp:
    delta_r = 0.45e-7
    uH = 303

    tempNumerator = lambda t: np.exp(-2 * t / tauSlow)
    numerator = quad(tempNumerator, 0, tSlow)

    tempDenominator = lambda u: hx.specific_heat(u) / hx.resistivity_BG(u, rrr=rrr, rho273=r273, delta_r=delta_r)
    denominator = quad(tempDenominator, uH, uf)

    A = I0 * np.sqrt(numerator[0] / (rho * denominator[0]))

    print('With warm end temp at {} K and an increased resistivity of {:.1e}'.format(uH, delta_r))
    print('Minimum Burnout Area is {:.3f} cm^2'.format(A))
    print('Minimum Burnout Diameter is {:.3f} cm'.format(np.sqrt(4 * A / np.pi)))