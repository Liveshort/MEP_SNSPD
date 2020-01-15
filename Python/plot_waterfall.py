import numpy as np
import matplotlib.pyplot as plt
import struct

with open("../sim_results/param.info") as file:
    paramList = []
    for line in file:
        paramList.extend(line.split("; ")[1:])

    runtype = int(paramList.pop(0))
    J0 = int(paramList.pop(0))
    J1 = int(paramList.pop(0))
    N = int(paramList.pop(0))
    timeskip = int(paramList.pop(0))
    ETratio = int(paramList.pop(0))
    numberOfT = int(paramList.pop(0))
    numberOfI = int(paramList.pop(0))
    numberOfR = int(paramList.pop(0))
    numberOfC = int(paramList.pop(0))
    dX0 = float(paramList.pop(0))
    dX1 = float(paramList.pop(0))
    dt = float(paramList.pop(0))

print(dX0, dX1, dt)

T0 = []
T1 = []
I1 = []
I2 = []
R0 = []
R1 = []
V_c = []

with open("../sim_results/T.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*J0*N//timeskip)):
        T0.append(item)
    for (item, ) in struct.iter_unpack('d', file.read(8*J1*N//timeskip)):
        T1.append(item)

T0 = np.array(T0).reshape(N//timeskip, J0)
T1 = np.array(T1).reshape(N//timeskip, J1)

with open("../sim_results/I.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        I1.append(item)
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        I2.append(item)

I1 = np.array(I1)
I2 = np.array(I2)

with open("../sim_results/R.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        R0.append(item)
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        R1.append(item)

R0 = np.array(R0)
R1 = np.array(R1)

with open("../sim_results/V_c.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
        V_c.append(item)

V_c = np.array(V_c)

t = np.arange(0, N, timeskip)
tE = np.arange(0, N*ETratio, timeskip//10)
x0 = np.linspace(-(J0//2 - 1)*dX0, J0//2*dX0, J0)
x1 = np.linspace(-(J1//2 - 1)*dX1, J1//2*dX1, J1)

plt.figure()
plt.contourf(t*dt*1e9, x0*1e6, T0.transpose() + 1e-12, 50, cmap="hot")
plt.xlabel("t (ns)")
plt.ylabel("x ($\mu$m)")
plt.colorbar(label="T (K)")
plt.show(block=False)

plt.figure()
plt.contourf(t*dt*1e9, x1*1e6, T1.transpose() + 1e-12, 50, cmap="hot")
plt.xlabel("t (ns)")
plt.ylabel("x ($\mu$m)")
plt.colorbar(label="T (K)")
plt.show(block=False)

plt.figure()
plt.plot(tE*dt*1e9, I1*1e6)
plt.xlabel("t (ns)")
plt.ylabel("I1 ($\mu$A)")
plt.show(block=False)

plt.figure()
plt.plot(tE*dt*1e9, I2*1e6)
plt.xlabel("t (ns)")
plt.ylabel("I2 ($\mu$A)")
plt.show(block=False)

plt.figure()
plt.plot(tE*dt*1e9, R0/1000)
plt.xlabel("t (ns)")
plt.ylabel("R0 (wire) (k$\Omega$)")
plt.show(block=False)

plt.figure()
plt.plot(tE*dt*1e9, R1/1000)
plt.xlabel("t (ns)")
plt.ylabel("R1 (wire) (k$\Omega$)")
plt.show()
