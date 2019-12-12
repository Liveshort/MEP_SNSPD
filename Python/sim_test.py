import numpy as np
import ctypes as ct

so_file = "../C/snspd.so"
snspd = ct.CDLL(so_file)

class SimRes(ct.Structure):
    _fields_ = [
        ("J",               ct.c_size_t),
        ("N",               ct.c_size_t),
        ("T",               ct.POINTER(ct.c_char)),
        ("I",               ct.POINTER(ct.c_char)),
        ("R",               ct.POINTER(ct.c_char)),
        ("exitValue",       ct.c_int)
    ]

class SimData(ct.Structure):
    _fields_ = [
        # general info
        ("J",               ct.c_size_t),
        ("N",               ct.c_size_t),
        ("numberofI",       ct.c_size_t),
        ("numberofR",       ct.c_size_t),
        # physical dimensions
        ("wireLength",      ct.c_double),
        ("wireThickness",   ct.c_double),
        ("wireWidth",       ct.c_double),
        ("tMax",            ct.c_double),
        # experiment specific data
        ("T_c",             ct.c_double),
        ("I_c0",            ct.c_double),
        ("c_p",             ct.c_double),
        ("c_e",             ct.c_double),
        ("alpha",           ct.c_double),
        ("T_sub",           ct.c_double),
        # data specific to the standard model (runtype 0)
        ("R_L_std",         ct.c_double),
        ("C_m_std",         ct.c_double),
        ("I_b_std",         ct.c_double),
        ("initHS_l_std",    ct.c_double),
        ("initHS_T_std",    ct.c_double),
        ("rho_norm_std",    ct.c_double),
        ("L_w_std",         ct.c_double)
    ]

simData = SimData()
simData.J = ct.c_size_t(20)
simData.N = ct.c_size_t(10)
simData.numberOfI = ct.c_size_t(2)
simData.numberOfR = ct.c_size_t(2)
simData.wireLength = ct.c_double(5.6e-6)

print(simData)

res = snspd.run_snspd_simulation(simData, ct.c_int(0))

print(res)
