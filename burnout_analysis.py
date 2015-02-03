__author__ = 'Chris Jones'
__email__ = 'cjones@tantaya.com'

"""
For a cryogenic current carrying lead, there exists a minimum area that will provide protection from over heating. Determining the amount of temperature change in a given time, under given conditions is the goal of the following analyses.

Our given conditions include some material specific constants, some expected project values, and conditions we want to simulate.

The basic question that this analysis answers is: How much area does a conductor have to have so that it won't go above a certain temperature?

We ask this question under a variety of circumstances. For example, with an exponentially decaying conductor current. Or, what if the conductor is warmer than usual to begin with? Or, what effect would cold working the conductor, which increases its resistivity, have on the minimum area?

The material constants are specific to C10200 OFHC copper, as are the functions used from heat_xfer.py. 

Units used:
gram
centimeter
second

Copyright (c) 2015 Christopher W Jones

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import numpy as np  # numpy v1.9.1
from scipy.integrate import quad  # scipy v0.14.0
import matplotlib.pyplot as plt  # matplotlib v1.4.2
import heat_xfer as hx

# Given Analysis Conditions
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
r_c_constant = 1
slow_dump = 1
fast_dump = 1
cw = 1
temp_var = 1
cw_and_temp = 1

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

    print('In a slow, {:.1f} s, dump:'.format(tSlow))
    print('With starting temp {} K, and ending temp {} K'.format(u0, uf))
    print('\nUsing temp dependent specific heat and resistivity:\nMinimum area: {:.4f} cm^2'.format(A))
    print('Minimum diameter: {:.4f} cm'.format(np.sqrt(4 * A / np.pi)))

# Using temp dependent specific heat and resistivity for fast dump
if fast_dump:
    tempNumerator = lambda t: np.exp(-2 * t / tauFast)
    numerator = quad(tempNumerator, 0, tFast)

    tempDenominator = lambda u: hx.specific_heat(u) / hx.resistivity_BG(u, rrr=rrr, rho273=r273)
    denominator = quad(tempDenominator, u0, uf)

    A = I0 * np.sqrt(numerator[0] / (rho * denominator[0]))

    print('\n\nIn a fast, {:.1f} s, dump:'.format(tFast))
    print('With starting temp {} K, and ending temp {} K'.format(u0, uf))
    print('\nUsing temp dependent specific heat and resistivity:\nMinimum area: {:.4f} cm^2'.format(A))
    print('Minimum diameter: {:.4f} cm'.format(np.sqrt(4 * A / np.pi)))

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
    # print(A)
    b = np.sqrt(4 * A / np.pi)
    # print(b)
    # print(A[-1]-A[0])
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

    # print(A)

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