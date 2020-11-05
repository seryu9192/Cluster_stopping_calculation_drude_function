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
def drude_function(x):
    #parameters from Z.Tan et al, Rad.Env.Biophys. 45 2(2006) 135-143
    Z_ave = 4 #for Glycine
    a = 0.2882925381282788 # defined to satisfy f-sum rule
    b = (19.927 + 0.9807 * Z_ave) / E_0
    c = (13.741 + 0.3215 * Z_ave) / E_0
    return a*x/((x**2-b**2)**2 + (c*x)**2)
