from math import pi, sqrt, cos, sin
from jinja2 import Undefined
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

M = 1
P = 20
X = 40
X_S = 2
Y = 40
T = 50
D = 1000


# ODE equations adapted from jasmcole at https://jasmcole.com/2014/09/04/black-holes-and-fractal-basins/
def odeBH(state, t, m):
    x = state[0]
    y = state[1]
    ux = state[2]
    uy = state[3]
    xbh = 0
    ybh = 0
    M = m

    if sqrt(x**2 + y**2) <= 2*M:
        return [0, 0, 0, 0]

    U = BlackHolePotential(xbh, ybh, M, x, y)
    [dUdx, dUdy] = GradientBlackHolePotential(xbh, ybh, M, x, y)
    udotgradU = ux*dUdx + uy*dUdy

    dxdt = ux/U
    dydt = uy/U
    duxdt = ((1 + 2*(ux**2 + uy**2))*dUdx - ux*udotgradU)/(U**2)
    duydt = ((1 + 2*(ux**2 + uy**2))*dUdy - uy*udotgradU)/(U**2)

    return [dxdt, dydt, duxdt, duydt]


def BlackHolePotential(xbh, ybh, M, x, y):
    u = 1
    return u + M/sqrt((x - xbh)**2 + (y - ybh)**2)


def GradientBlackHolePotential(xbh, ybh, M, x, y):
    dUdx = 0
    dUdy = 0

    dUdx = dUdx - M*(x - xbh)*((x - xbh) ** 2 + (y - ybh)**2)**-1.5
    dUdy = dUdy - M*(y - ybh)*((x - xbh) ** 2 + (y - ybh)**2)**-1.5
    return [dUdx, dUdy]


def draw(event=Undefined):
    t = np.linspace(0, T, D)
    trajectories = []
    for i in range(P):
        #u0 = [X_S*i+round(-P/2), 0, 10, -1]
        u0 = [-20, 0, cos(0.03*i), sin(0.03*i)]
        trajectories.append(
            odeint(odeBH, u0, t, (M,), mxstep=5000, rtol=1.4e-1, atol=1.4e-5, hmin=0.1))

    plt.xlim(-X/2, X/2)
    plt.ylim(-Y/2, Y/2)

    plt.grid()
    # solarRing = plt.Circle((0, 0), 3*M, color="grey", alpha=0.3)
    # eventHorizon = plt.Circle((0, 0), 2*M, color="black")
    ax.plot(0, 0, marker="o", markersize=1,
            markerfacecolor="black", markeredgecolor="black")
    # ax.add_patch(solarRing)
    # ax.add_patch(eventHorizon)
    for p in trajectories:
        ax.plot(p[:, 0], p[:, 1], color='orange')


def remove():
    for i, line in enumerate(ax.lines):
        line.remove()
        ax.lines.pop(i)
        # ax.lines[i].remove()
    ax.lines[0].remove()
    ax.lines[0].remove()


def photonUpdate(x):
    global P
    P = x
    draw(Undefined)


fig, ax = plt.subplots()
fig.set_size_inches(6, 6)
# fig.subplots_adjust(bottom=0.4)

draw()

plt.ylabel('y')
plt.xlabel('x')

plt.show()
