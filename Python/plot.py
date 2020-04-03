import numpy as np
import matplotlib.pyplot as plt
import struct

from math import ceil

# returns number of currents for a given runtype
def rtToNC(runtype):
    if runtype == 0:
        return 4
    elif runtype == 1:
        return 5
    elif runtype == 4:
        return 7
    elif runtype == 6:
        return 9

# returns number of bias currents (or thermal branches, if you want) for a given runtype
def rtToNBC(runtype):
    if runtype == 0:
        return 1
    elif runtype == 1:
        return 1
    elif runtype == 4:
        return 2
    elif runtype == 6:
        return 3

with open("../sim_results/param.info") as file:
    paramList = []
    for line in file:
        paramList.extend(line.split("; ")[1:])

    runtype = int(paramList.pop(0))
    trnstype = int(paramList.pop(0))

    J = []
    for i in range(rtToNBC(runtype)):
        J.append(int(paramList.pop(0)))

    N = int(paramList.pop(0))
    timeskip = int(paramList.pop(0))
    ETratio = int(paramList.pop(0))
    numberOfT = int(paramList.pop(0))
    numberOfI = int(paramList.pop(0))
    numberOfR = int(paramList.pop(0))
    numberOfC = int(paramList.pop(0))

    I_b = list()
    dX = list()

    for i in range(numberOfT):
        I_b.append(float(paramList.pop(0)))
    for i in range(numberOfT):
        dX.append(float(paramList.pop(0)))

    dt = float(paramList.pop(0))

print("dX: {} [m], dt: {} [s]".format(dX, dt))
#print(runtype, trnstype, J, N, timeskip, ETratio, numberOfT, numberOfI, numberOfR, numberOfC, I_b, dX, dt)

T_tmp = [[] for i in range(numberOfT)]
I = [[] for i in range(numberOfI)]
R = [[] for i in range(numberOfR)]
V_c = [[] for i in range(numberOfC)]

with open("../sim_results/T.bin", "rb") as file:
    for i in range(numberOfT):
        for (item, ) in struct.iter_unpack('d', file.read(8*J[i]*ceil(N/timeskip))):
            T_tmp[i].append(item)

T0 = np.array(T_tmp[0]).reshape(ceil(N/timeskip), J[0])
if numberOfT > 1:
    T1 = np.array(T_tmp[1]).reshape(ceil(N/timeskip), J[1])
if numberOfT > 2:
    T2 = np.array(T_tmp[2]).reshape(ceil(N/timeskip), J[2])

with open("../sim_results/I.bin", "rb") as file:
    for i in range(numberOfI):
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            I[i].append(item)

I = np.array(I)

with open("../sim_results/R.bin", "rb") as file:
    for i in range(numberOfR):
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            R[i].append(item)

R = np.array(R)

with open("../sim_results/V_c.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        V_c.append(item)

V_c = np.array(V_c)

t = np.arange(0, N-1, timeskip)
tE = np.arange(0, N*ETratio-1, timeskip//10)

x0 = np.linspace(-(J[0]//2 - 1)*dX[0], J[0]//2*dX[0], J[0])
if numberOfT > 1:
    x1 = np.linspace(-(J[1]//2 - 1)*dX[1], J[1]//2*dX[1], J[1])
if numberOfT > 2:
    x2 = np.linspace(-(J[2]//2 - 1)*dX[2], J[2]//2*dX[2], J[2])

# plot the thermal results
vmin, vmax = (2.7, 12.0)
if numberOfT > 2:
    fig, axes = plt.subplots(nrows=3, ncols=1)
    im = axes.flat[0].contourf(t*dt*1e9, x0*1e6, T0.transpose() + 1e-12, 50, cmap="hot", vmin=vmin, vmax=vmax)
    axes.flat[0].set_ylabel("x ($\mu$m)")
    axes.flat[1].contourf(t*dt*1e9, x1*1e6, T1.transpose() + 1e-12, 50, cmap="hot", vmin=vmin, vmax=vmax)
    axes.flat[1].set_ylabel("x ($\mu$m)")
    axes.flat[2].contourf(t*dt*1e9, x2*1e6, T2.transpose() + 1e-12, 50, cmap="hot", vmin=vmin, vmax=vmax)
elif numberOfT > 1:
    fig, axes = plt.subplots(nrows=2, ncols=1)
    im = axes.flat[0].contourf(t*dt*1e9, x0*1e6, T0.transpose() + 1e-12, 50, cmap="hot", vmin=vmin, vmax=vmax)
    axes.flat[0].set_ylabel("x ($\mu$m)")
    axes.flat[1].contourf(t*dt*1e9, x1*1e6, T1.transpose() + 1e-12, 50, cmap="hot", vmin=vmin, vmax=vmax)
else:
    fig, axes = plt.subplots(nrows=1, ncols=1)
    im = axes.flat[0].contourf(t*dt*1e9, x0*1e6, T0.transpose() + 1e-12, 50, cmap="hot", vmin=vmin, vmax=vmax)
plt.xlabel("t (ns)")
plt.ylabel("x ($\mu$m)")

im.set_clim(vmin, vmax)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
plt.ylabel("T [K]")

plt.show(block=False)

# plot the currents
plt.figure()
for i in range(numberOfI - 2):
    plt.plot(tE*dt*1e9, I[i,:]*1e6)
if runtype == 0:
    plt.legend(["I0 (detector)", "I_load"], loc="upper center")
elif runtype == 1:
    plt.legend(["I0 (detector)", "I1 (par stg 0)", "I_load"], loc="upper center")
elif runtype == 4:
    plt.legend(["I0 (detector)", "I1 (par stg 0)", "I2 (wire stg 1)", "I3 (par stg 1)", "I_load"], loc="upper center")
elif runtype == 6:
    plt.legend(["I0 (detector)", "I1 (par stg 0)", "I2 (wire stg 1)", "I3 (par stg 1)", "I4 (wire stg 2)", "I5 (par stg 2)", "I_load"], loc="upper center")
plt.xlabel("t (ns)")
plt.ylabel("I ($\mu$A)")
plt.show(block=False)

# plot the currents in the transmission line
plt.figure()
for i in range(numberOfI - 3, numberOfI):
    plt.plot(tE*dt*1e9, I[i,:]*1e6)
plt.legend(["I_load", "I_tl_sum", "I_tl"], loc="upper center")
plt.xlabel("t (ns)")
plt.ylabel("I ($\mu$A)")
plt.show(block=False)

plt.figure()
plt.plot(tE*dt*1e9, R[0]/1000)
if numberOfR > 1:
    plt.plot(tE*dt*1e9, R[1]/1000)
if numberOfR > 2:
    plt.plot(tE*dt*1e9, R[2]/1000)
plt.xlabel("t (ns)")
plt.ylabel("R (wire) (k$\Omega$)")
plt.show()
