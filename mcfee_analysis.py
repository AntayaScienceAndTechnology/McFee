__author__ = 'Chris Jones'

import numpy as np
import McFee_Functions as mf
import matplotlib.pyplot as plt

# Given conditions
uL = 40.0
uH = 293.0
du = 0.1
num = (uH - uL + du) * (1 / du)
num = int(round(num))
rrr = 150
rho273 = 1.71e-6
I = 200.

# Graphs available (Set desired flag to '1')
current_deviation = 0  # Graph showing how deviations from design current affect heat load
optimized_la_ratio = 0  # Graph of l to a ratio for given conditions
optimized_heat_load = 0  # Graph of optimizied heat load by design current
cold_work = 0  # Simulate added resistivity due to cold work done at 295 K
temp_variation = 0  # Simulation variation in temp of lead at warm end
cw_and_temp = 0  # Combine the cold work and temp variation simulation

if current_deviation:
    x = np.logspace(-3, 3, 50, base=2)
    q_qmin = 0.5 * (x + 1 / x)
    plt.semilogx(x, q_qmin, 'k')
    plt.ylim([0, 4])
    plt.xscale('log', basex=2)
    plt.grid(which='minor', color='k')
    plt.grid(which='major', color='k')
    plt.xlabel('Ratio of Actual Current to Design Current')
    plt.ylabel(r'Ratio of $Q$ to $Q_{min}$')
    plt.title('Relationship Between Minimum Heat Load and Design Current')
    plt.show()

if optimized_la_ratio:
    u = np.linspace(uH, uL - 10, num=num)
    calc_du = u[0] - u[1]
    temp_of_interest = uL  # This temp will be annotated on the graph
    toi_u = np.linspace(uH, uL, num=num)
    toi_du = toi_u[0] - toi_u[1]
    rat_interest = mf.get_la_ratio(toi_u, toi_du, I)
    print('Ratio from {} K to {} K is {:.1f}'.format(uH, uL, rat_interest))
    rat = np.zeros(len(u))

    for i in range(2, len(u) + 1):
        rat[i - 1] = mf.get_la_ratio(u[:i], calc_du, I)

    fig = plt.figure(num=0, figsize=(14, 11))
    ax = fig.gca()
    plt.ylim([0, 250])
    ax.annotate(r'Ratio at {} K = {:.1f}'.format(temp_of_interest, rat_interest),
                arrowprops=dict(facecolor='black', shrink=0.05, width=2),
                xytext=(100, 205), textcoords='data', size=20,
                xy=(temp_of_interest, rat_interest), xycoords='data',
                bbox=dict(boxstyle='square', fc='white'))
    plt.plot(u, rat, 'k')
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

    plt.plot(I_range, y, 'k')
    plt.title('Optimized Current Lead Heat Load')
    plt.ylabel('Heat Load per Lead (W)')
    plt.xlabel('Lead Current (A)')
    plt.grid()
    plt.show()

if cold_work:
    r25 = 0.3e-7
    r50 = 0.45e-7
    u = np.arange(uH, uL - du, -du)
    delta_r = np.linspace(0., 0.6e-7, num=25)

    ratio = np.zeros(len(delta_r))
    rat25 = mf.get_la_ratio(u, du, I, 0.3e-7)
    rat50 = mf.get_la_ratio(u, du, I, 0.45e-7)

    for i in range(len(delta_r)):
        ratio[i] = mf.get_la_ratio(u,du, I, delta_r[i])

    fig = plt.figure()
    ax = fig.gca()
    ax.annotate('25% CW\n{:.1f}'.format(rat25),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2),
                    xytext=(1e-8, rat25 - 1), textcoords='data', size=20,
                    xy=(r25, rat25), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    ax.annotate('50% CW\n{:.1f}'.format(rat50),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2),
                    xytext=(2.5e-8, rat50 - 1), textcoords='data', size=20,
                    xy=(r50, rat50), xycoords='data',
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

    heat25 = mf.get_qn2(uL, uH, I, du, 0.3e-7)
    heat50 = mf.get_qn2(uL, uH, I, du, 0.45e-7)

    fig = plt.figure()
    ax = fig.gca()
    ax.annotate('25% CW\n{:.2f}'.format(heat25),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2),
                    xytext=(0.3e-7, 8.82), textcoords='data', size=20,
                    xy=(r25, heat25), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    ax.annotate('50% CW\n{:.2f}'.format(heat50),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2),
                    xytext=(4.5e-8, 8.9), textcoords='data', size=20,
                    xy=(r50, heat50), xycoords='data',
                    bbox=dict(boxstyle='square', fc='white'))
    plt.plot(delta_r, heat, 'k')
    plt.grid(True)
    plt.xlabel(r'Increase in Resistivity $(\Omega cm)$')
    plt.ylabel(r'Minimum Heat Load at Cold End (W)')
    plt.title(r'Cold Working Effect on Minimum Heat Load at Cold End')
    plt.show()

if temp_variation:
    k = 16
    delta_u = np.linspace(280., 305., num=k)
    heat = np.zeros(k)
    ratio = np.zeros(k)

    for i in range(len(delta_u)):
        da = np.linspace(delta_u[i],uL, num)
        du = da[0] - da[1]
        ratio[i] = mf.get_la_ratio(da, du, I)
        heat[i] = mf.get_qn2(uL, delta_u[i], I, du)

    plt.figure()
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
