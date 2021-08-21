# common_library.py

from numpy import pi

#CONSTANTS(atomic unit)
a_0 = 0.529 #Bohr radius(A)
E_0 = 27.2 #Hartree energy(eV)
e2 = 14.4 #(eV A)
E_PROTON = 25
Z_PROTON = 1 #atomic number of proton
E_HELIUM = 100 #(keV at 1 au velocity)
Z_HELIUM = 2 #atomic number of helium
E_CARBON = 300 #(keV at 1 au velocity)
Z_CARBON = 6 #atomic number of carbon
A = 0.240 # constant for calculating screening length of Brandt-Kitagawa model

#properties of solid amino acid
AMINO_PROP = {
    "Gly":{
        "Z_ave": 4.0,
        "a": 0.2882925381282788,
        "Ep": 19.6487041570916, 
        "vF": 1.069795654
    },
    "Phe":{
        "Z_ave": 3.8260869565217392,
        "a": 0.2882925381282788,
        "Ep": 20.4041164494357,
        "vF": 1.097042407
    },
}

#constants for Thomas-Fermi-Moliere screening function
TMF = {
    "alpha": [0.10, 0.55, 0.35],
    "beta": [6.0, 1.20, 0.30]
}

#convert index to index string (i, j) => "ij" (i < j) or "ji" (i > j) 
def str_ind(i, j):
    res = str(i) + str(j)
    if res > res[::-1]:res = res[::-1]
    return res

#Drude-type Optical energy loss function
def drude_function(x, amino_name):
    #parameters from Z.Tan et al, Rad.Env.Biophys. 45 2(2006) 135-143
    Z_ave = AMINO_PROP[amino_name]["Z_ave"]
    a = AMINO_PROP[amino_name]["a"] # defined to satisfy f-sum rule
    b = (19.927 + 0.9807 * Z_ave) / E_0 #to atomic unit
    c = (13.741 + 0.3215 * Z_ave) / E_0 #to atomic unit
    return a*x/((x**2-b**2)**2 + (c*x)**2)

# Energy loss function for
def elf_low_v(k, w, amino_name):
    E_p = AMINO_PROP[amino_name]["Ep"]/E_0
    v_F = AMINO_PROP[amino_name]["vF"]
    k_TF = 3**(1/2)*E_p/v_F # Thomas-Fermi
    return pi/(2*v_F)*(k_TF**2)*k*w/((k**2+k_TF**2)**2)
