#!/usr/bin/env python3

import numpy as np
from scipy import special
import sys
import scipy as sp
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import ROOT
from scipy.interpolate import PchipInterpolator

mV2V   = 1e-3
um2cm  = 1e-4
MeV2eV = 1e6
W_eV   = 19.5
alpha0 = 0.21
keV2eV   = 1000.0
param_set = {"Scalettar": [3.9, 3.0], "Aprile": [3.7, 2.8], 'Ereditato': [4.2, 2.3]}

element_charge_C  = 1.602176634e-19

#Ekin: keV, F: V/cm
def r_Segreto(Ekin, F, param):
    Wi     = 23.6 #eV
    k     = param[0] * mV2V  #V
    gamma = param[1] / um2cm #/cm

    alpha = 0.227 * (MeV2eV**2)         #eV^2/cm
    beta  = 1.7   *  MeV2eV             #eV/cm
    Ee    = 1 + k * gamma / F
    
    z  = Ekin * keV2eV * (F + k * beta / Wi) / (k * alpha / Wi)
    R  = ( F / (F + k * beta / Wi) ) * (1 - np.log(1 + z) / z)
    Re = R * Ee
    
    return Re

# #from C to keV
# def convert_Q2keV(Q, F, param):
#     E_array = np.linspace(1,10000, 999)
#     r = r_Segreto(E_array, F, param)

#     Q_array = convert_keV2Q(E_array, r)
#     func_Q2E = sp.interpolate.interp1d(Q_array, E_array, kind='cubic')
#     return func_Q2E(Q)

#from C to keV
def convert_Q2keV_tatsumi(Q, F, param):
    E_array = np.linspace(1,10000, 10000)
    r = r_Segreto(E_array, F, param)

    Q_array = convert_keV2Q(E_array, r)

    interp = PchipInterpolator(Q_array, E_array, extrapolate=False)

    E = interp(Q)
    Qmin, Qmax = Q_array[0], Q_array[-1]
    Q = np.asarray(Q, float)
    E = np.where(~np.isfinite(E) & (Q <= Qmin), E_array[0], E)
    E = np.where(~np.isfinite(E) & (Q >= Qmax), E_array[-1], E)
    return E

    

def convert_keV2Q(E_keV, r):
    return E_keV * keV2eV * r / (W_eV * (1 + alpha0)) * element_charge_C

def hill(E, A, n, Ec, p):
    return A * E**n / (1.0 + (Ec/E)**p)

def sbpl(E, A, a1, a2, Eb, Delta):
    x = E / Eb
    return A * x**a1 * (1.0 + x**(1.0/Delta))**((a2 - a1)*Delta)

def log_parabola(Q, a, b, c):
    return np.exp(a + b*np.log(Q) + c*np.log(Q)**2)

def sbpl_inv(Q, Eb, beta1, beta2, Qb, Delta):
    x = Q / Qb
    return Eb * x**beta1 * (1.0 + x**(1.0/Delta))**((beta2-beta1)*Delta)


if __name__=="__main__":
    E_array = np.logspace(0, 4, 99)
    F_array = np.array([200, 300, 500, 1000, 2000, 3000])
    model_recombination = "Aprile"

    numFieldValue = len(F_array)
    
    Q_array = np.zeros((F_array.shape[0], E_array.shape[0]))
    
    tgraphList = []
    spl3List = []
    
    for i, F in enumerate(F_array):
        r = r_Segreto(E_array, F, param_set[model_recombination])
        Q_array[i] = convert_keV2Q(E_array, r)
        #spline補完
        tgraphList.append(ROOT.TGraph(len(E_array), Q_array[i], E_array))
        spl3List.append(ROOT.TSpline3(f"spl{i}", tgraphList[i]))
        spl3List[i].SetName(f"E{F}")
    
    parList = []
    f = ROOT.TFile("QvsEkeV_spline.root", "recreate")
    
    for i in range(numFieldValue):
        p0 = (3e2, 0.55, 0.95, 1e-15, 0.3)
        pars, cov = curve_fit(sbpl_inv, Q_array[i], E_array, p0=p0, maxfev=20000)
        spl3List[i].Write()
        parList.append(pars)
    
    f.Close()
