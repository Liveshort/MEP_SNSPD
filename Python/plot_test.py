import numpy as np
import matplotlib.pyplot as plt
import struct

with open("../sim_results/param.info") as file:
    paramList = []
    for line in file:
        paramList.append(line.split("; ")[1])

    J = int(paramList[0])
    N = int(paramList[1])
    timeskip = int(paramList[2])
    numberOfT = int(paramList[3])
    numberOfI = int(paramList[4])
    numberOfR = int(paramList[5])
    dX = float(paramList[6])
    dt = float(paramList[7])

print(dX, dt)

T = []
I = []
R = []

with open("../sim_results/T.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*J*N//timeskip)):
        T.append(item)

T = np.array(T).reshape(N//timeskip, J)

with open("../sim_results/I.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N//timeskip)):
        I.append(item)

I = np.array(I)

with open("../sim_results/R.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N//timeskip)):
        R.append(item)

R = np.array(R)

t = np.arange(0, N-1, timeskip)
x = np.linspace(-(J//2 - 1)*dX, J//2*dX, J)

plt.contourf(t*dt*1e9, x*1e6, T.transpose() + 1e-12, 50, cmap="hot")
plt.xlabel("t (ns)")
plt.ylabel("x ($\mu$m)")
plt.colorbar()
plt.show(block=False)

plt.figure()
plt.plot(t*dt*1e9, I*1e6)
plt.xlabel("t (ns)")
plt.ylabel("I ($\mu$A)")
plt.show(block=False)

plt.figure()
plt.plot(t*dt*1e9, R/1000)
plt.xlabel("t (ns)")
plt.ylabel("R (wire) (k$\Omega$)")
plt.show()
