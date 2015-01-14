__author__ = 'Chris Jones'

import numpy as np
import McFee_Functions as mf
import matplotlib.pyplot as plt

# Given conditions
uL = 40
uH = 293
du = 1
rrr = 150
rho273 = 1.71e-6
I = 200.

# Graphs available
current_deviation = 0  # Graph showing how deviations from design current affect heat load
optimized_la_ratio = 0  # Graph of l to a ratio for given conditions
optimized_heat_load = 0  # Graph of optimizied heat load by design current
cold_work = 0
temp_variation = 0
cw_and_temp = 1

if current_deviation:
    x = np.logspace(-3, 3, 50, base=2)
    q_qmin = 0.5 * (x + 1 / x)
    plt.semilogx(x, q_qmin)
    plt.ylim([0, 4])
    plt.xscale('log', basex=2)
    plt.grid(which='minor', color='k')
    plt.grid(which='major', color='k')
    plt.xlabel('Ratio of Actual Current to Design Current')
    plt.ylabel(r'Ratio of $Q$ to $Q_{min}$')
    plt.title('Relationship Between Minimum Heat Load and Design Current')
    plt.show()

if optimized_la_ratio:
    u = np.arange(uH, uL - du, -du)
    cell_temp = mf.get_cell_temps(u)
    rat40 = mf.get_la_ratio(np.arange(uH, 40 - du, -du), du, I)
    rat = np.zeros(len(u))

    for i in range(2, len(u)):
        rat[i] = mf.get_la_ratio(u[0:i], du, I)

    fig = plt.figure(num=0, figsize=(14, 11))
    ax = fig.gca()
    ax.annotate(r'Ratio at 40 K = {:.1f}'.format(rat40),
                arrowprops=dict(facecolor='black', shrink=0.02, width=2),
                xytext=(100, 205), textcoords='data', size=20,
                xy=(40, rat40), xycoords='data',
                bbox=dict(boxstyle='square', fc='white'))
    plt.plot(cell_temp, rat)
    plt.title('Optimized Current Lead Length to Area Ratio')
    plt.ylabel('Ratio (L/A)')
    plt.xlabel('Temperature (K)')
    plt.grid()
    plt.show()


if optimized_heat_load:
    I_range = np.arange(150, 250, 0.5)
    y = np.zeros(len(I_range))

    for i in range(len(I_range)):
        y[i] = mf.get_qn2(uL, uH, I_range[i], du)

    plt.plot(I_range, y)
    plt.title('Optimized Current Lead Heat Load')
    plt.ylabel('Heat Load per Lead (W)')
    plt.xlabel('Lead Current (A)')
    plt.grid()
    plt.show()

if cold_work:
    u = np.arange(uH, uL - du, -du)
    delta_r = np.linspace(0., 0.6e-7, num=25)

    ratio = np.zeros(len(delta_r))
    rat25 = mf.get_la_ratio(u, du, I, 0.3e-7)
    rat50 = mf.get_la_ratio(u, du, I, 0.45e-7)

    for i in range(len(delta_r)):
        ratio[i] = mf.get_la_ratio(u,du, I, delta_r[i])

    print(ratio)

    fig = plt.figure()
    ax = fig.gca()
    ax.annotate('25% CW\n{:.1f}'.format(rat25),
                    arrowprops=dict(facecolor='black', shrink=0.02, width=2),
                    xytext=(1e-8, 207), textcoords='data', size=20,
                    xy=(0.3e-7, 207.7), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    ax.annotate('50% CW\n{:.1f}'.format(rat50),
                    arrowprops=dict(facecolor='black', shrink=0.02, width=2),
                    xytext=(4e-8, 208.2), textcoords='data', size=20,
                    xy=(0.45e-7, 206.7), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    plt.plot(delta_r, ratio, 'k')
    plt.grid(True)
    plt.xlabel(r'Increase in Resistivity $(\Omega cm)$')
    plt.ylabel(r'Optimized Conductor (L/A) Ratio')
    plt.title(r'Cold Working Effect on Conductor Length to Area Ratio')
    plt.show()

    heat = np.zeros(len(delta_r))
    for i in range(len(delta_r)):
        heat[i] = mf.get_qn2(uL, uH, I, du, delta_r[i])

    print(heat)

    heat25 = mf.get_qn2(uL, uH, I, du, 0.3e-7)
    heat50 = mf.get_qn2(uL, uH, I, du, 0.45e-7)

    fig = plt.figure()
    ax = fig.gca()
    ax.annotate('25% CW\n{:.2f}'.format(heat25),
                    arrowprops=dict(facecolor='black', shrink=0.02, width=2),
                    xytext=(0.1e-7, 9), textcoords='data', size=20,
                    xy=(0.29e-7, 8.94), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    ax.annotate('50% CW\n{:.2f}'.format(heat50),
                    arrowprops=dict(facecolor='black', shrink=0.02, width=2),
                    xytext=(4e-8, 8.9), textcoords='data', size=20,
                    xy=(0.45e-7, 9.02), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    plt.plot(delta_r, heat, 'k')
    plt.grid(True)
    plt.xlabel(r'Increase in Resistivity $(\Omega cm)$')
    plt.ylabel(r'Minimum Heat Load at Cold End (W)')
    plt.title(r'Cold Working Effect on Minimum Heat Load at Cold End')
    plt.show()

if temp_variation:
    k = 50
    delta_u = np.linspace(280, 305, num=k)
    du = 1e-1
    heat = np.zeros(k)
    ratio = np.zeros(k)



    for i in range(len(delta_u)):
        da = np.arange(delta_u[i], 40., -du)
        ratio[i] = mf.get_la_ratio(da, du, I)
        heat[i] = mf.get_qn2(40., delta_u[i], I, du)

    plt.plot(delta_u, ratio, 'k')
    plt.grid(True)
    plt.xlabel(r'Temperature of Warm End of Conductor (K)')
    plt.ylabel(r'Optimized Conductor (L/A) Ratio')
    plt.title('Temperature Variation at Warm End of Conductor\'s\nEffect on Conductor Length to Area Ratio')
    plt.show()

    plt.plot(delta_u, heat, 'k')
    plt.grid(True)
    plt.xlabel(r'Temperature of Warm End of Conductor (K)')
    plt.ylabel(r'Minimum Heat Load at Cold End (W)')
    plt.title('Temperature Variation at Warm End of Conductor\'s\nEffect on Minimum Heat Load at Cold End')
    plt.show()

if cw_and_temp:
    delta_r = 0.45e-7
    uH = 303
    u = np.arange(uH, uL - du, -du)

    ratio = mf.get_la_ratio(u, du, I, delta_r)
    heat = mf.get_qn2(uL, uH, I, du, delta_r)

    print('With warm end temp at {} K and an increased resistivity of {:.1e}'.format(uH, delta_r))
    print('Minimum Heat Load at Cold End is {:.2f}\nOptimum Length to Area Ratio is {:.1f}'.format(heat, ratio))
