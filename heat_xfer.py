__author__ = 'Chris Jones'
__maintainer__ = 'Chris Jones'
__email__ = 'cjones@tantaya.com'
__status__ = 'Development'

"""
heat_xfer.py

Contains functions for temperature dependent physical parameters of OFHC copper.
Each function has its source listed in the docstring. There are valid ranges of temperature for each function.

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
import numpy as np  #v1.9.1
from scipy.integrate import quad  #v0.14.0

DENSITY_CU = 8.96  # g / cm^3
INTR_RESIST_CU = 1.71e-6  # Ohm * cm
ALPHA_CU = 3.86e-3  # 1/K ... Should this be 6.8e-3


def therm_cond_cu(t, rrr=150):
    """
        Returns thermal conductivity of OFHC copper for
        a given t in Kelvin and purity ratio rrr.
        Returned values have units of W/(cm K).
        Valid between 4 K and 300 K only.
        t must be an integer or numpy array
    :rtype : float
        Valid for the following rrr values: 50, 100, 150
        Source http://cryogenics.nist.gov/MPropsMAY/OFHC%20Copper/OFHC_Copper_rev.htm
    """
    singleVal = 0
    if not hasattr(t, '__iter__'):
        t = np.arange(t, t + 1, 1)
        singleVal = 1

    a = {50: 1.8743, 100: 2.2154, 150: 2.3797}
    b = {50: -0.41538, 100: -0.47461, 150: -0.4918}
    c = {50: -0.6018, 100: -0.88068, 150: -0.98615}
    d = {50: 0.13294, 100: 0.13871, 150: 0.13942}
    e = {50: 0.26426, 100: 0.29505, 150: 0.30475}
    f = {50: -0.0219, 100: -0.02043, 150: -0.019713}
    g = {50: -0.051276, 100: -0.04831, 150: -0.046897}
    h = {50: 0.0014871, 100: 0.001281, 150: 0.0011969}
    i = {50: 0.003723, 100: 0.003207, 150: 0.0029988}

    output = (1 / 100) * 10 ** ((a[rrr] + c[rrr] * t ** 0.5 + e[rrr] * t + g[rrr] * t ** 1.5 + i[rrr] * t ** 2) /
                                (1 + b[rrr] * t ** 0.5 + d[rrr] * t + f[rrr] * t ** 1.5 + h[rrr] * t ** 2))

    if singleVal:
        return output[0]

    return output


def resistivity_BG(t, rrr=150, rho273=1.71e-6, a=9.323e-6, delta_r=0.):
    """
    Adapted from:
    Theoretical and Mathematical Physics, 166(1): 37–42 (2011)
    THE BLOCH–GRUNEISEN FUNCTION OF ARBITRARY ORDER AND ITS SERIES REPRESENTATIONS
    Valid 10 mK to 1357 K
    a, rho, and rrr are valid for OFHC C10200 copper.
    Return units: Ohm cm
    :param t:
    :param rrr:
    :param rho273:
    :return:
    """

    scalar_value = 0
    if not hasattr(t, '__iter__'):
        t = np.arange(t, t + 1, 1)
        scalar_value = 1

    rho_o = rho273 / rrr  # Ohm cm

    debyeT = 343.5  # K
    n = 5

    output = np.zeros(len(t))

    intg = lambda x: (x ** n) / ((np.exp(x) - 1) * (1 - np.exp(-x)))

    for el in range(len(t)):
        T = t[el]
        x, y = quad(intg, 0, debyeT / T,)
        output[el] = (T / debyeT) ** n * x

    if scalar_value:
        return rho_o + delta_r + a * output[0]

    output *= a

    output += rho_o + delta_r

    return output  # Ohm cm


def specific_heat(b):
    """
    Taken from:
    http://www.nist.gov/data/PDFfiles/jpcrd263.pdf
    pp. 1253-1254
    :param b:
    :return:
    """
    singleVal = 0

    if not hasattr(b, '__iter__'):
        b = np.arange(b, b + 1, 1)
        singleVal = 1

    a = []
    c = []
    temp = b / 100.0

    for u in temp:

        if u < 0.25:
            raise ValueError('Temperature must be above 25 K')
        elif 0.25 <= u < 0.29:
            a.extend([0.96297, 12.4259, 39.9045, 64.7626])
            u_min = 0.25
        elif 0.29 <= u < 0.4473:
            a.extend([1.528, 15.929, 47.676, -89.7484])
            u_min = 0.29
        elif u < 0.692:
            a.extend([4.864, 24.266, 5.3237, -29.051])
            u_min = 0.4473
        elif u < 1.3946:
            a.extend([10.695, 21.653, -16.003, 5.2425])
            u_min = 0.692
        elif u < 2.0:
            a.extend([19.827, 6.93, -4.9524, 1.8736])
            u_min = 1.3946
        elif u < 3.3:
            a.extend([22.623, 2.9936, -1.5496, 0.3724])
            u_min = 2.0
        elif u < 12.37:
            a.extend([24.714, 0.8526, -0.09737, 0.00873])
            u_min = 3.3
        else:
            raise ValueError('Temperature must be below 1237 Kelvin.')

        u -= u_min

        c_J_per_mol_K = a[0] + a[1] * u + a[2] * u ** 2 + a[3] * u ** 3

        c_J_per_gram_K = c_J_per_mol_K / 63.546  # Atomic Weight Cu = 63.546 g/mol

        c.extend([c_J_per_gram_K])
        a = []

    if singleVal:
            return c[0]

    return c


def specific_heat1(t):
    """
        Returns specific heat of OFHC copper for a given t in Kelvin.
        Returned values have units of J/(g*K).
        Valid between 4 K and 300 K only.
        Source http://cryogenics.nist.gov/MPropsMAY/OFHC%20Copper/OFHC_Copper_rev.htm
    """
    # Polynomial Coefficients
    a = -1.91844
    b = -0.15973
    c = 8.61013
    d = -18.996
    e = 21.9661
    f = -12.7328
    g = 3.54322
    h = -0.3797

    # Polynomial
    return (1 / 1000) * 10 ** (a + b * np.log10(t) + c * np.log10(t) ** 2 + d * np.log10(t) ** 3
                         + e * np.log10(t) ** 4 + f * np.log10(t) ** 5 + g * np.log10(t) ** 6 + h * np.log10(t) ** 7)


def therm_cond_al(t):
    """
        Returns thermal conductivity of 6061-T6 Aluminum for
        a given t in Kelvin.
        Returned values have units of W/(cm K).
        Valid between 4 K and 300 K only.
        t must be an integer or numpy array
    :rtype : float
        Source http://cryogenics.nist.gov/Papers/Cryo_Materials.pdf
    """
    a = 0.07918
    b = 1.09570
    c = -0.07277
    d = 0.08084
    e = 0.02803
    f = -0.09464
    g = 0.04179
    h = -0.00571

    # if type(t) == np.ndarray and (t.min() < 4 or t.max() > 300):
    #     return 0
    #
    # if (type(t) == int or type(t) == float) and (t < 4 or t > 300):
    #     return 0

    return (1 / 100) * 10 ** (a + b * np.log10(t) + c * np.log10(t) ** 2 + d * np.log10(t) ** 3
                              + e * np.log10(t) ** 4 + f * np.log10(t) ** 5 + g * np.log10(t) ** 6 + h * np.log10(t) ** 7)


def expansion_coeff(t):
    """
        Returns expansion coefficient of OFHC copper for a given t in Kelvin.
        Returned values have units of 1/K.
        Valid between 4 K and 300 K only.
        Source http://cryogenics.nist.gov/MPropsMAY/OFHC%20Copper/OFHC_Copper_rev.htm
    """
    # Polynomial Coefficients
    a = -17.9081289
    b = 67.131914
    c = -118.809316
    d = 109.9845997
    e = -53.8696089
    f = 13.30247491
    g = -1.30843441

    # if not 4 <= t.all() <= 300:
    #     return 0

    # Polynomial
    return (1 / 1000) * 10 ** (a + b * np.log10(t) + c * np.log10(t) ** 2 + d * np.log10(t) ** 3
                               + e * np.log10(t) ** 4 + f * np.log10(t) ** 5 + g * np.log10(t) ** 6)