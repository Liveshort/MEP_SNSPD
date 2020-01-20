# SNSPD (Superconducting Nanowire Single-Photon Detector) Simulation Code
## Running the simulation (Linux)
* First clone the repository to a place of your liking:
```bash
git clone git@github.com:Liveshort/MEP_SNSPD.git
```
* Install some required packages (if I missed one, you'll get a `package not found` error somewhere down the line, just install that package analogous to the ones here):
```bash
sudo apt install build-essential gfortran python3-dev python3-pip
```
* Install required Python3 packages:
```bash
pip3 install numpy matplotlib
```
* Download the current LAPACK version (version 3.9.0 was used at the time of writing) from http://www.netlib.org/lapack/.
* Unpack the LAPACK repository by opening a terminal, going to the downloads folder and running the following commands:
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
```
* Four `lib[***].a` files will have appeared in your LAPACK folder, copy those to the empty C/lib folder in the copy of the git repo on your pc
* Open a terminal, move into the C folder and run the following command to run a simulation:
```bash
make all && time make run args="../sim_setup/setup_yang.info ../sim_results/"
```
* If all went well, everything will compile and a progress bar will appear. It should be done in a few seconds for a simple simulation. The first argument is the input for the simulation, found in the sim_setup folder, and the second argument is the output folder, in this case the sim_results/ folder. Some `[***].bin` files will have appeared here, along with an info file containing information about the simulation in plain text (which you can read as well).
* To look at the results of your simulation in nice figures, open the `Python/plot_[***].py` file in your favorite editor/IDE and run the code there or with the following command:
```bash
python3 plot_[***].py
```
* To run other simulations provided by this simulation software, provide different setup files with the `args` argument in the `make run` command.
## Running the simulation (Windows)
* Running the simulations under Windows is slightly more complicated than under Linux, but is the same for the most part. It does require the Windows Subsystem for Linux (WSL).
* Install the Windows Subsystem for Linux (https://docs.microsoft.com/en-us/windows/wsl/install-win10), choose any distribution you like, but know that the code was tested under Ubuntu-like distributions, so 16.04 or the newer 18.04 are the safest choices.
* Open a Windows Command Prompt and enter `bash`. This will open a terminal in WSL.
* Now follow the steps of the Linux preparations above, except the Python code part. Note that the `ulimit` command might fail, which caused the LAPACK test suite to fail on my machine. The libraries, however, compiled just fine, so I could just copy them over to their designated folder. This will probably also be the case for you.
* You can run Python natively on your Windows machine. Open the Python folder in the copy of this repo on your pc in PyCharm, Spyder or another Python interpreter of your liking.
* Hit `run`, or the equivalent in your software.
