import numpy as np
import matplotlib.pyplot as plt
import struct

with open("../simres/param.info") as file:
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

with open("../simres/T.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*J*N//timeskip)):
        T.append(item)

T = np.array(T).reshape(N//timeskip, J)

with open("../simres/I.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N//timeskip)):
        I.append(item)

I = np.array(I)

with open("../simres/R.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N//timeskip)):
        R.append(item)

R = np.array(R)

plt.figure()
plt.imshow(T)
plt.colorbar()
plt.show(block=False)

plt.figure()
plt.plot(I)
plt.ylabel("I (A)")
plt.show(block=False)

plt.figure()
plt.plot(R)
plt.ylabel("R (wire) (Ohm)")
plt.show()
