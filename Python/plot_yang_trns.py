import numpy as np
import matplotlib.pyplot as plt
import struct

with open("../sim_results/param.info") as file:
    paramList = []
    for line in file:
        paramList.append(line.split("; ")[1])

    runtype = int(paramList.pop(0))
    trns = int(paramList.pop(0))
    J = int(paramList.pop(0))
    N = int(paramList.pop(0))
    timeskip = int(paramList.pop(0))
    ETratio = int(paramList.pop(0))
    numberOfT = int(paramList.pop(0))
    numberOfI = int(paramList.pop(0))
    numberOfR = int(paramList.pop(0))
    numberOfC = int(paramList.pop(0))
    I_b = float(paramList.pop(0))
    dX = float(paramList.pop(0))
    dt = float(paramList.pop(0))

print(dX, dt)

T = []
I = [[] for i in range(4+runtype)]
R = []
V_c = []

with open("../sim_results/T.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*J*N//timeskip)):
        T.append(item)

T = np.array(T).reshape(N//timeskip, J)

with open("../sim_results/I.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        I[0].append(item)
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        I[1].append(item)
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        I[2].append(item)
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        I[3].append(item)
    if runtype == 1:
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            I[4].append(item)

I = np.array(I)

with open("../sim_results/R.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        R.append(item)

R = np.array(R)

with open("../sim_results/V_c.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        V_c.append(item)

V_c = np.array(V_c)

t = np.arange(0, N-1, timeskip)
tE = np.arange(0, N*ETratio-1, timeskip//10)
x = np.linspace(-(J//2 - 1)*dX, J//2*dX, J)

plt.contourf(t*dt*1e9, x*1e6, T.transpose() + 1e-12, 50, cmap="hot")
plt.xlabel("t (ns)")
plt.ylabel("x ($\mu$m)")
plt.colorbar(label="T (K)")
plt.show(block=False)

plt.figure()
plt.plot(tE*dt*1e9, I[0,:]*1e6)
plt.plot(tE*dt*1e9, I[1,:]*1e6)
plt.plot(tE*dt*1e9, I[2,:]*1e6)
plt.plot(tE*dt*1e9, I[3,:]*1e6)
if runtype == 1:
    plt.plot(tE*dt*1e9, I[4,:]*1e6)
plt.xlabel("t (ns)")
plt.ylabel("I ($\mu$A)")
plt.legend(["I0 (detector)", "Iload", "Itlsum", "Itl"], loc="upper center")
plt.ylim([0, 20])
plt.show(block=False)

plt.figure()
plt.plot(tE*dt*1e9, R/1000)
plt.xlabel("t (ns)")
plt.ylabel("R (wire) (k$\Omega$)")
plt.show()
