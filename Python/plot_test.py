import numpy as np
import matplotlib.pyplot as plt
import struct

J = 1000
N = 10000

T = []
I = []
R = []

with open("../simres/T.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*J*N)):
        T.append(item)

T = np.array(T).reshape(N, J)

with open("../simres/I.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N)):
        I.append(item)

I = np.array(I)

with open("../simres/R.bin", "rb") as file:
    for (item, ) in struct.iter_unpack('d', file.read(8*N)):
        R.append(item)

R = np.array(R)

plt.figure()
plt.imshow(T.transpose())
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
