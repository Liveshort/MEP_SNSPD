# SNSPD (Superconducting Nanowire Single-Photon Detector) Simulation Code
## Setting up the simulation (Linux)
* First install some required packages (if I missed one, you'll get a `package not found` error somewhere down the line, just install that package analogous to the ones here):
```bash
sudo apt install git build-essential gfortran python3-dev python3-pip
```
* Then clone the repository to a place of your liking:
```bash
git clone git@github.com:Liveshort/MEP_SNSPD.git
```
* Install required Python3 packages:
```bash
pip3 install numpy matplotlib
```
* Download either the current OpenBLAS library from https://github.com/xianyi/OpenBLAS (recommeded, because it is much more optimized and runs multicore) [Follow the steps below 1.] **OR** the current LAPACK version (version 3.9.0 was used at the time of writing) from http://www.netlib.org/lapack/ [Follow the steps below 2.].
### ...using OpenBLAS (recommended)
* Clone the git repository of OpenBLAS somewhere on your system by opening a terminal, `cd`-ing into a folder of your liking and running the following commands:
```bash
# clone repository
git clone https://github.com/xianyi/OpenBLAS.git
# move into folder
cd OpenBLAS
# make a folder to put the output files
mkdir OUT
# make BLAS, LAPACK, LAPACKE
make PREFIX=OUT/ install
```
* A `libopenblas_***-r***.dev.a` file will have appeared in the OUT/lib folder. Copy this file to the C/lib folder in the SNSPD project repository and rename it to `libopenblas.a`. Header files should already be included, but if need be, you can copy over your own generated files from the OUT/include folder into the C/include_openblas folder. Note that the header files have been changed slightly (all references to complex.h have been commented out) due to a conflict of the variable name `I`. This should not bother you if you just use the supplied `.h` files.
### ...using Netlib LAPACK
* Unpack the LAPACK repository by opening a terminal, going to the Downloads folder and running the following commands:
```bash
# unpack
tar -xvf lapack-3.9.0.tar.gz
# move into folder and copy config file
cd lapack-3.9.0
cp make.inc.example make.inc
# increase the stack size during compilation of the LAPACK.
# it will probably compile without this command, but the test suite will not run.
ulimit -s unlimited
# make BLAS and LAPACK
make all
# move into LAPACKE folder and make that too
cd LAPACKE
make
# move into CBLAS folder and make that too
cd ../CBLAS
make
```
* Five `lib[***].a` files will have appeared in your LAPACK folder, copy those to the empty C/lib folder in the copy of the git repo on your PC.
## Running the simulation (Linux)
* Open a terminal, move into the C folder and run the following command to run a simulation (OpenBLAS or Netlib will be selected automatically):
```bash
make all && time make run args="../sim_setup/setup_yang.info ../sim_results/"
```
* If all went well, everything will compile and a progress bar will appear. It should be done in a few seconds for a simple simulation. The first argument is the input for the simulation, found in the sim_setup folder, and the second argument is the output folder, in this case the sim_results/ folder. Some `[***].bin` files will have appeared here, along with an info file containing information about the simulation in plain text (which you can read as well).
* To look at the results of your simulation in nice figures, open the `Python/plot_[***].py` file in your favorite editor/IDE and run the code there or with the following command:
```bash
python3 plot_[***].py
```
* To run other simulations provided by this simulation software, provide different setup files with the `args` argument in the `make run` command. Understanding what happens is easiest by looking in the Docs folder for supporting images of the electrical circuits used. To understand the data model, it is best to look at the setup files, and the Python scripts supplied here. Results are saved in binary format for maximum speed, where the order of the data can be understood by taking a look at the source code of the Python scripts.
## Setting up and running the simulation (Windows, using Windows Subsystem for Linux)
* Running the simulations under Windows is slightly more complicated than under Linux, but is the same for the most part. It does require the Windows Subsystem for Linux (WSL).
* Install the Windows Subsystem for Linux (https://docs.microsoft.com/en-us/windows/wsl/install-win10), choose any distribution you like, but know that the code was tested under Ubuntu-like distributions, so 16.04 or the newer 18.04 are the safest choices.
* Open a Windows Command Prompt and enter `bash`. This will open a terminal in WSL.
* Now follow the steps of the Linux preparations above, except the Python code part. Note that, when using the Netlib LAPACK implementation, the `ulimit` command might fail, which caused the LAPACK test suite to fail on my machine. The libraries, however, compiled just fine, so I could just copy them over to their designated folder. This will probably also be the case for you.
* You can run Python natively on your Windows machine. Open the Python folder in the copy of this repo on your pc in PyCharm, Spyder or another Python interpreter of your liking.
* Hit `run`, or the equivalent in your software.
