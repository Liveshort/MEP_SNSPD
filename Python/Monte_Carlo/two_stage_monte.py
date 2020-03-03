import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import struct

from subprocess import PIPE

def write_setup(name, values, emptyLines, intLines):
    with open(name, 'w+') as fp:
        currLine = 1
        for val in values:
            if currLine in emptyLines:
                fp.write("    // skip\n")
                currLine += 1
            if currLine in intLines:
                fp.write("{}; $\n".format(int(val)))
            else:
                fp.write("{}; $\n".format(val))
            currLine += 1

def fill_values(input):
    valuesList = []

    for i in range(1, input.shape[1]):
        valuesList.append(list(input[str(i)]))

    return valuesList

def read_results(outputFolder, I, R, it):
    global runtype, J0, J1, N, timeskip, ETratio, numberOfT, numberOfI, numberOfR, numberOfC, I_b0, I_b1, dX0, dX1, dt

    for i in range(4):
        I[i].append([])
    for i in range(2):
        R[i].append([])

    with open("{}/param.info".format(outputFolder)) as file:
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
        I_b0 = float(paramList.pop(0))
        I_b1 = float(paramList.pop(0))
        dX0 = float(paramList.pop(0))
        dX1 = float(paramList.pop(0))
        dt = float(paramList.pop(0))

    with open("{}/I.bin".format(outputFolder), "rb") as file:
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            I[0][it].append(item)
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            I[1][it].append(item)
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            I[2][it].append(item)
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            I[3][it].append(item)

    with open("{}/R.bin".format(outputFolder), "rb") as file:
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            R[0][it].append(item)
        for (item, ) in struct.iter_unpack('d', file.read(8*N*ETratio//timeskip*10)):
            R[1][it].append(item)

    return

if __name__ == "__main__":
    # indicate lines that should be empty and lines that should be integers
    twoStageEmptyLines = [3, 10, 14, 18, 25, 45]
    twoStageIntLines = [1, 2, 4, 5, 6, 8, 9]

    # read input from csv
    input = pd.read_csv("setup.csv")

    # make a list of values
    valuesList = fill_values(input)

    # set up the target I and R arrays
    I = list()
    R = list()

    for i in range(4):
        I.append([])
    for i in range(2):
        R.append([])

    # define some globals from the sims that come in handy
    global runtype, J0, J1, N, timeskip, ETratio, numberOfT, numberOfI, numberOfR, numberOfC, I_b0, I_b1, dX0, dX1, dt

    # run simulations
    for it, values in enumerate(valuesList):
        write_setup("setup.info", values, twoStageEmptyLines, twoStageIntLines)

        print(subprocess.run(["../../C/bld/simulation", "setup.info", "res/"]))

        read_results("res/", I, R, it)

    # plot simulations
    tE = np.arange(0, N*ETratio, timeskip//10)

    # extract currents
    I0 = np.transpose(np.array(I[0]))
    I1 = np.transpose(np.array(I[1]))
    I2 = np.transpose(np.array(I[2]))
    I3 = np.transpose(np.array(I[3]))
    R0 = np.transpose(np.array(R[0]))
    R1 = np.transpose(np.array(R[1]))

    plt.figure(figsize=(11,6))
    plt.subplot(2,3,1)
    for i in range(5):
        plt.plot(tE*dt*1e9, I0[:,i]*1e6)
    plt.xlabel("t (ns)")
    plt.ylabel("I ($\mu$A)")
    plt.title("I$_0$ (nanowire current)")
    plt.subplot(2,3,2)
    for i in range(5):
        plt.plot(tE*dt*1e9, I2[:,i]*1e6)
    plt.xlabel("t (ns)")
    plt.ylabel("I ($\mu$A)")
    plt.title("I$_2$ (stg 1 current)")
    plt.subplot(2,3,3)
    for i in range(5):
        plt.plot(tE*dt*1e9, I3[:,i]*1e6)
    plt.xlabel("t (ns)")
    plt.ylabel("I ($\mu$A)")
    plt.title("I$_3$ (stg 1 par current)")
    plt.subplot(2,3,4)
    for i in range(5):
        plt.plot(tE*dt*1e9, (I_b0 + I_b1 - I0[:,i]- I1[:,i]- I2[:,i]- I3[:,i])*1e6)
    plt.xlabel("t (ns)")
    plt.ylabel("I ($\mu$A)")
    plt.title("I$_{res}$ (meas current)")
    plt.subplot(2,3,5)
    for i in range(5):
        plt.plot(tE[tE*dt*1e9<1]*dt*1e9, R0[tE*dt*1e9<1,i]/1e3)
    plt.xlabel("t (ns)")
    plt.ylabel("R (k$\Omega$)")
    plt.title("R$_0$ (nanowire res)")
    plt.subplot(2,3,6)
    for i in range(5):
        plt.plot(tE[tE*dt*1e9<1]*dt*1e9, R1[tE*dt*1e9<1,i]/1e3)
    plt.xlabel("t (ns)")
    plt.ylabel("R (k$\Omega$)")
    plt.title("R$_1$ (stg 1 res)")
    plt.tight_layout(pad=1, w_pad=0.9, h_pad=0.8)
    plt.show()
