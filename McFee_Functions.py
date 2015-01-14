__author__ = 'Chris Jones'
import sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, 'C:/Users/Chris Jones/Documents/PythonProjects/OFHC_Therm_Cond')
import heat_xfer as hx


def get_kr_ave1(u_l, u_h):
    """
    Returns the average k*r value from u_l to u_h using integration
    :param u_l:
    :param u_h:
    :return:
    """
    if u_l == u_h:
        return 0.

    kr_func = lambda u: hx.therm_cond_cu(u, rrr=150) * hx.resistivity_BG(u)
    return (1 / (u_h - u_l)) * quad(kr_func, u_l, u_h)[0]


def get_qn1(u_l, u_h, I):  # The heat flow at the cold end
    """
    Returns the nth heat flow value in watts using integration
    :param u_l:
    :param u_h:
    :param I:
    :return:
    """
    return I * np.sqrt(2 * get_kr_ave1(u_l, u_h) * (u_h - u_l))


def get_qps1(u, u_l,  I=200):
    qps = np.zeros(len(u))

    for i in range(len(u)):
        qps[i] = get_qn1(u_l, u[i], I)

    return qps


def get_kr_ave2(u_l, u_h, du, r_increase=0.):
    """
    Returns the average k*r value from u_l to u_h using summation
    :param u_l:
    :param u_h:
    :param du:
    :return:
    """
    if u_l == u_h:
        return 0.

    x = np.arange(u_l, u_h + du, du)

    return np.sum(get_kp(x)*get_rp(x, r_increase=r_increase)) / len(x)


def get_qn2(u_l, u_h, I, du, r_increase=0.):
    """
    Returns the nth heat flow value in watts using summation
    :param u_l: cold end temperature in Kelvin
    :type u_l: float
    :param u_h: warm end temperature in Kelvin
    :type u_h: float
    :param I: conductor current
    :type I: float
    :param du: temperature step size
    :type du: float
    :returns: heat load seen at cold end of conductor
    :rtype: float
    """
    return I * np.sqrt(2 * get_kr_ave2(u_l, u_h, du, r_increase=r_increase) * (u_h - u_l))


def get_qps2(u, du, I=200):
    qps = np.zeros(len(u))
    u_max = max(u)

    for i in range(len(u)):
        qps[i] = get_qn2(u_max, u[i], I, du)

    return qps


def get_kr_cumsum(cell_temps, r_increase=0.):
    """
    For a given cell temperature range, return the sum of the k*r products.
    Used by 'get_la_ratio()'
    :param cell_temps:
    :type cell_temps: numpy.ndarray
    :return: cumulative sum of k*r values
    :rtype: numpy.ndarray
    """
    return np.cumsum(get_kp(cell_temps)*get_rp(cell_temps, r_increase=r_increase))


def get_qn3(kr_sum, du, I):
    """
    Returns the value of Q from a range. Meant to be used by 'get_la_ratio()'
    :param kr_sum:
    :param du:
    :param I:
    :return:
    :rtype:
    """
    return I * np.sqrt(2 * kr_sum * du)


def get_la_ratio(u, du, I, r_increase=0.):
    """
    Given a temperature range and current, returns the optimized length to area ratio of the conductor.
    :param r_increase:
    :param u:
    :param du:
    :param I:
    :return:
    """
    ct = get_cell_temps(u)
    kr = get_kr_cumsum(ct, r_increase=r_increase)
    sp = get_sp(ct, r_increase=r_increase)

    qs = get_qn3(kr, du, I)

    ratio = 0.

    for i in range(len(ct) - 1):
        ratio += (sp[i] - sp[i+1]) * qs[i+1]

    ratio += sp[-1] * qs[-1]

    return ratio / I ** 2


def get_kp(u, rrr=150):
    """
    Given a temperature, or array of temperatures, returns the thermal conductivity.
    :param u:
    :param rrr:
    :return:
    """
    return hx.therm_cond_cu(u, rrr)


def get_sp(u, rrr=150, rho273=1.71e-6, r_increase=0.):
    """
    Given a temperature, or an array of temperatures, returns the electrical conductivity.
    :param r_increase:
    :param u:
    :param rrr:
    :param rho273:
    :return:
    """
    return 1 / (hx.resistivity_BG(u, rrr, rho273) + r_increase)


def get_rp(u, rrr=150, rho273=1.71e-6, r_increase=0.):
    """
    Given a temperature, or an array of temperatures, returns the electrical resistivity.
    :param u:
    :param rrr:
    :param rho273:
    :return:
    """
    return hx.resistivity_BG(u, rrr, rho273) + r_increase


def get_cell_temps(u):
    """
    For a temperature array of length n >= 2, where the values represent temperatures at cell boundaries,
    this function returns the average temperature for the cell in an n-1 length array.
    :param u:
    :return:
    """
    temps = np.zeros(len(u)-1)

    for i in range(len(temps)):
        temps[i] = 0.5 * (u[i] + u[i + 1])

    return temps