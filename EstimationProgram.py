""" Sebastian Strug
Last updated: 26/06/2017
Estimating thickness of 2D materials based on optical microscope images 

Largely based on the work presented in:
"Making graphene visible"
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
from scipy import interpolate
from pylab import meshgrid,cm,imshow,colorbar
import sys
import os
from tkFileDialog import askopenfilename
import time
import glob
import re

"""We will be looking at the visible rangle of light, 400nm-780nm"""
lda = np.linspace(0.4,0.78,1000) #lda stands for lambda, create 1000 evenly spaced values

"""
Extract the experimental data for refractive indices (real and imaginary),
at different wavelengths for: Silicon (pure), hexagonal Boron Nitride, Antimony
""" 
scriptDir = sys.path[0] #this specifies the directory we are in
Si_handle = open(os.path.join(scriptDir, "Refractive indices/refracIndex_Si.txt"), 'r') #open a handle to the file, ospath joins the directory with another path
Sidata = [[float(x) for x in line.split()] for line in Si_handle.readlines()] # Do a double-nested list comprehension to get the data from read in lines into a matrix
ldaSi = [Sidata[i][0] for i in range(len(Sidata))] #first column of the data, wavelength
Sinreal = [Sidata[i][1] for i in range(len(Sidata))] #second column, real refractive index
Sinimag = [Sidata[i][2] for i in range(len(Sidata))] #third column of the data, imaginary refractive index
Si_2_handle = open(os.path.join(scriptDir, "Refractive indices/refracIndex_Si_2.txt"), 'r') 
Sidata_2 = [[float(x) for x in line.split()] for line in Si_2_handle.readlines()[1:]] #need to skip the header in this file
ldaSi_2 = [Sidata_2[i][2] for i in range(len(Sidata_2))]
Sinreal_2 = [Sidata_2[i][3] for i in range(len(Sidata_2))]
Sinimag_2 = [Sidata_2[i][4] for i in range(len(Sidata_2))]
SiO2_handle = open(os.path.join(scriptDir, "Refractive indices/refracIndex_SiO2.txt"), 'r')
SiO2data = [[float(x) for x in line.split()] for line in SiO2_handle.readlines()]
ldaSiO2 = [SiO2data[i][0] for i in range(len(SiO2data))] #for silicon dioxide
SiO2nreal = [SiO2data[i][1] for i in range(len(SiO2data))]
patternForData = re.compile("^\s+|\s*,\s*|\s+$") #this is for all the whitespaces and commas in the data (data was scraped using 'Data Thief' on a graph from a paper)
BN_handle = open(os.path.join(scriptDir, "Refractive indices/refracIndex_hBN.txt"), 'r')
BNdata = [[float(x) for x in patternForData.split(line) if x] for line in BN_handle.readlines()] #slightly different since this isn't tab deliminated
ldaBN = [BNdata[i][0] for i in range(len(BNdata))] #for hexagonal Boron Nitride
BNnreal = [BNdata[i][1] for i in range(len(BNdata))]
BNnimag = [BNdata[i][2] for i in range(len(BNdata))]
BP_real_handle = open(os.path.join(scriptDir, "Refractive indices/refracIndexReal_BP"), 'r')
BP_real_data = [[float(x) for x in patternForData.split(line) if x] for line in BP_real_handle.readlines()[1:]] #another 'Data Thief' program scrape from a paper
ldaBPreal = [BP_real_data[i][0]/1000 for i in range(len(BP_real_data))] #for hexagonal Boron Nitride
BPnreal = [BP_real_data[i][1] for i in range(len(BP_real_data))]
BP_imag_handle = open(os.path.join(scriptDir, "Refractive indices/refracIndexImag_BP"), 'r')
BP_imag_data = [[float(x) for x in patternForData.split(line) if x] for line in BP_imag_handle.readlines()[1:]] #another 'Data Thief' program scrape from a paper
ldaBPimag = [BP_imag_data[i][0]/1000 for i in range(len(BP_imag_data))] #for hexagonal Boron Nitride
BPnimag = [BP_imag_data[i][1] for i in range(len(BP_imag_data))]
Sb_handle = open(os.path.join(scriptDir, "Refractive indices/refracIndex_Sb.txt"), 'r')
Sbdata = [[float(x) for x in patternForData.split(line) if x] for line in Sb_handle.readlines()]
ldaSb = [Sbdata[i][0] for i in range(len(Sbdata))] #for Antimony
Sbnreal = [Sbdata[i][1] for i in range(len(Sbdata))]
Sbnimag = [Sbdata[i][2] for i in range(len(Sbdata))]
WSe2Path = os.path.join(scriptDir, "Refractive indices/WSe2_data.csv")
WSe2data = np.genfromtxt(WSe2Path, delimiter=',')
ldaWSe2 = WSe2data[0]/1000
WSe2nreal = WSe2data[1]
WSe2nimag = WSe2data[2]
MoS2Path = os.path.join(scriptDir, "Refractive indices/MoS2_data.csv")
MoS2data = np.genfromtxt(MoS2Path, delimiter=',')
ldaMoS2 = MoS2data[0]/1000
MoS2nreal = MoS2data[1]
MoS2nimag = MoS2data[2]
NbSe2Path = os.path.join(scriptDir, "Refractive indices/NbSe2_data.csv")
NbSe2data = np.genfromtxt(NbSe2Path, delimiter=',')
ldaNbSe2 = NbSe2data[0]/1000
NbSe2nreal = NbSe2data[1]
NbSe2nimag = NbSe2data[2]

"""Interpolate these values so we have refractive indices for any wavelength in
the visible light range we want"""
Sirealn = interpolate.interp1d(ldaSi, Sinreal) #real refractive index, n
Siimagn = interpolate.interp1d(ldaSi, Sinimag) #imaginary refractive index, k
Sirealn_2 = interpolate.interp1d(ldaSi_2, Sinreal_2) #real refractive index, n
Siimagn_2 = interpolate.interp1d(ldaSi_2, Sinimag_2) #imaginary refractive index, k
SiO2realn = interpolate.interp1d(ldaSiO2, SiO2nreal)
BNrealn = interpolate.interp1d(ldaBN, BNnreal)
BNimagn = interpolate.interp1d(ldaBN, BNnimag)
BPrealn = interpolate.interp1d(ldaBPreal, BPnreal)
BPimagn = interpolate.interp1d(ldaBPimag, BPnimag)
Sbrealn = interpolate.interp1d(ldaSb, Sbnreal)
Sbimagn = interpolate.interp1d(ldaSb, Sbnimag)  
WSe2realn = interpolate.interp1d(ldaWSe2, WSe2nreal)
WSe2imagn = interpolate.interp1d(ldaWSe2, WSe2nimag)
MoS2realn = interpolate.interp1d(ldaMoS2, MoS2nreal)
MoS2imagn = interpolate.interp1d(ldaMoS2, MoS2nimag)
NbSe2realn = interpolate.interp1d(ldaNbSe2, NbSe2nreal)
NbSe2imagn = interpolate.interp1d(ldaNbSe2, NbSe2nimag)

"""
We calculate the refractive indices of Silicon Oxide and Sodium Chloride from fitted
Sellmeier's coefficients. These fits have been found to agree very well with experimental data.
"""

def refractiveindexNaCl(lda):
    #returns n for NaCl, in visible light range it is transparent so no k (extinguishing factor), only n
    #from H. H. Li, Journal of Physical Chemistry Reference data, 5, 329, 1976
    #"Refractive index of alkali halides and its wavelength and temperature derivatives"
    A = 1.00055; B = 0.19800; C = 0.48398
    D = 0.38696; E = 0.25998; F = 0.08796
    G = 3.17064; H = 0.30038
    nsquared = A + (B*lda**2)/(lda**2-0.05**2) + (C*lda**2)/(lda**2-0.1**2) + \
                 (D*lda**2)/(lda**2-0.128**2) + (E*lda**2)/(lda**2-0.158**2) + \
                 (F*lda**2)/(lda**2-40.5**2) + (G*lda**2)/(lda**2-60.98**2) + \
                 (H*lda**2)/(lda**2 - 120.34**2)
    n = np.sqrt(nsquared)
    return n      

"""For hBN, I found a paper that does not use the refractive index values of hBN
from the optical handbook (even though they cite that for SiO2/Si). Instead they use another
source and mention their data points which I replicate here. The paper is:
Optical thickness determination of hexagonal Boron Nitride flakes:"""
def hBNrefracLineVals(lda):
    tempx2 = float(560./1000)
    tempx1 = float(480./1000)
    tempy2 = float(1.85)
    tempy1 = tempy2*1.03
    tempgradient = (tempy2-tempy1)/(tempx2-tempx1)
    tempintercept = tempy2 - (tempgradient*tempx2)
    return ((tempgradient*lda) + tempintercept)

"""create 1000 valued refractive index array for each material, for each wavelength"""
#initialise our arrays
nvalsSi = np.zeros(len(lda),dtype=complex) 
nvalsSi_2 = np.zeros(len(lda),dtype=complex) 
nvalsBN = np.zeros(len(lda),dtype=complex)
nvalsBN_2 = np.zeros(len(lda),dtype=complex)
nvalsNaCl = np.zeros(len(lda),dtype=complex) 
nvalsGraphene = np.zeros(len(lda),dtype=complex)
nvalsGraphene_2 = np.zeros(len(lda),dtype=complex)
nvalsSb = np.zeros(len(lda),dtype=complex)
nvalsWSe2 = np.zeros(len(lda),dtype=complex)
nvalsMoS2 = np.zeros(len(lda),dtype=complex)
nvalsNbSe2 = np.zeros(len(lda),dtype=complex)
nvalsSiO2 = np.zeros(len(lda),dtype=complex)
nvalsMica = np.zeros(len(lda),dtype=complex)
nvalsMica_Muscovite = np.zeros(len(lda),dtype=complex)
nvalsMica_backtrack = np.zeros(len(lda),dtype=complex)
nvalsPDMS = np.zeros(len(lda),dtype=complex)
nvalsBP = np.zeros(len(lda),dtype=complex)

nvalsSi_test = np.zeros(len(lda),dtype=complex)
nvalsMica_test = np.zeros(len(lda),dtype=complex)

for i in range(0,len(lda)):
    nvalsSi[i] = Sirealn(lda[i]) + (Siimagn(lda[i])*(-1j)) #use interpolated data
    nvalsSi_2[i] = Sirealn_2(lda[i]) + (Siimagn_2(lda[i])*(-1j)) #use interpolated data
    nvalsSi = nvalsSi_2 #nvalsSi_2 is better
    nvalsSiO2[i] = SiO2realn(lda[i])-(0*1j) #from the imported SiO2 data 
    """Difference in SiO2 refractive indices is very minor"""
    nvalsBN[i] = BNrealn(lda[i]) - (BNimagn(lda[i])*(-1j)) #use interpolated data
    nvalsBN_2[i] = hBNrefracLineVals(lda[i])-(0*1j) #from the paper "Hunting for monolayer boron nitide"
    nvalsBN = nvalsBN_2 #change to the value in the paper
    nvalsNaCl[i] = refractiveindexNaCl(lda[i]) #from the function which uses Sellmeier coefficients
    nvalsGraphene[i] = 2.6-1.3j #refractive index of graphene is modelled as constant, from the paper
    nvalsGraphene_2[i] = 2.0 - 1.1j #"graphene thickness determination using reflection and contrast spectroscopy" paper says this is much better
    nvalsSb[i] = Sbrealn(lda[i]) + ( Sbimagn(lda[i])*(-1j) ) #use interpolated data
    nvalsMoS2[i] = MoS2realn(lda[i]) + ( MoS2imagn(lda[i])*(-1j) )
    nvalsWSe2[i] = WSe2realn(lda[i]) + ( WSe2imagn(lda[i])*(-1j) )
    nvalsNbSe2[i] = NbSe2realn(lda[i]) + ( NbSe2imagn(lda[i])*(-1j) )
    nvalsMica[i] = 2.51-(0*1j) #average of refractive index given on website
    #nvalsMica_Muscovite[i] = 1.55-(0*1j) #for muscovite mica from the paper "Atomically thin mica flakes and their application as ultrathin insulating substrates for graphene"
    nvalsMica_backtrack[i] = 1.19 #this is what the refractive index has to be roughly for the image detection contrast to work, for 542nm of light
    nvalsMica = nvalsMica_backtrack    
    nvalsPDMS[i] = 1.415 #from the paper "A new fabrication method for all-PDMS waveguides by Cai et al. (2013)
    #nvalsMica_test[i] = 1.82 #for 673nm wavelength of light
    if lda[i]<=min(ldaBPreal) or lda[i]<=min(ldaBPimag):
        nvalsBP[i] = BPrealn(ldaBPreal[0]) + (BPimagn(ldaBPimag[0])*(-1j))        
    else: nvalsBP[i] = BPrealn(lda[i]) + (BPimagn(lda[i])*(-1j))

nvals = np.zeros(len(lda),dtype=complex)
print("\n\nWhich material would you like to examine as your 2D material on top of a silicon substrate?\n\
If you want to enter your own refractive index, please type 'other': ")
materialType = raw_input("Fluorophlogopite, (h) Boron Nitride, Graphene, Antimony, Sodium Chloride, WSe2, MoS2, NbSe2, Black Phosphorus, or other? \n")
materialType = materialType.lower()
if materialType.lower() in ['mica', 'fluorophlogopite','m']:
    materialType = 'mica'
    nvals = nvalsMica
    formalName = 'Fluorophlogopite'
    monoLayer = 1./1000 #in um
elif materialType.lower() in ['hexagonal boron nitride','boron nitride','hbn','bn']:
    materialType = 'hbn'
    nvals = nvalsBN
    formalName = 'hexagonal Boron Nitride'
    monoLayer = 0.333/1000 #in nm
elif materialType.lower() in ['graphene','gr','g', 'grp']:
    materialType = 'grp'
    nvals = nvalsGraphene
    formalName = 'Graphene'
    monoLayer = 0.34/1000 #in nm
elif materialType.lower() in ['antimony','sb','a','the nicholas hemsworth exsberience']:
    nvals = nvalsSb    
    materialType = 'sb' 
    formalName = 'Antimony'
    monoLayer = 0.45/1000 #Thin NaCl films on silver (001): island growth and work function
elif materialType.lower() in ['sodium chloride', 'nacl', 'sc']:
    materialType = 'nacl'
    nvals = nvalsNaCl
    formalName = 'Sodium Chloride'
    monoLayer = 0.28/1000 #from paper "Thin NaCl films on silver (001): island growth and work function"
elif materialType.lower() in ['other','na']:
    refractiveIndexReal = input("Please enter the real component of the refractive index: ")
    refractiveIndexImag = input("Please enter the imaginary component of the refractive index: ")
    formalName = 'custom material'
    for i in range(0,len(lda)):
        nvals[i] = refractiveIndexReal + ( refractiveIndexImag*(-1j) )
elif materialType.lower() in ['wse2','wse','ws','w','ws2','tungsten','tungsten diselenide','tungstendiselenide']:
    nvals = nvalsWSe2    
    materialType = 'wse2' 
    formalName = 'WSe2'
    monoLayer = 0.67/1000 #Vahid says dichalcogenide monolayers are roughly 0.5nm
elif materialType.lower() in ['mos2','mos','ms','m','ms2','m2','molybdenum','molybdenum disulfide','molybdenumdisulfide']:
    nvals = nvalsMoS2    
    materialType = 'mos2' 
    formalName = 'MoS2'
    monoLayer = 0.675/1000
elif materialType.lower() in ['nbse2','nbse','nbs','nb','nb2','niobium diselenide','niobium','niobium selenide']:
    nvals = nvalsMoS2    
    materialType = 'nbse2' 
    formalName = 'NbSe2'
    monoLayer = 0.686/1000
elif materialType.lower() in ['bp','black phosphorus','black','phosphorus','blackphosphorus','bphos','bphosphorus']:
    nvals = nvalsBP
    materialType = 'bp'
    formalName = 'BP'
    monoLayer = 0.5/1000 #estimated, real value unknown
else:
    print('That is not a valid input, try again')    
    quit()    

print("Do you want to look at:\n \
1. Refractive index plot of {} \n \
2. Contrast as a function of wavelength for several SiO2 thicknesses for one {} thickness \n \
3. Contrast as a function of wavelength for several {} thicknesses for a SiO2 or PDMS substrate \n \
4. A colour plot of wavelength against SiO2 thickness with contrast as colour \n \
5. A colour plot of wavelength against {} thickness with contrast as colour \n \
6. Calculate the thickness for an image of {} I have the contrast for".format(formalName,formalName,formalName,formalName,formalName))
graphType = input("Please enter 1, 2, 3, 4, 5, or 6: ")
if graphType not in [1, 2, 3, 4, 5, 6, 'one','two','three','four', 'five','six']:
    print('That is not a valid input, try again')    
    quit
elif graphType in [2, 'two']:
    graphType = 2
    setThickness = float(input("Enter the thickness of the {} you would like to look at in nanometres: ".format(formalName)))
    setThickness = setThickness/1000   
    numberSiO2 = input("Enter the number of SiO2 thicknesses you would like to look at: ")
    SiO2List = []
    for i in xrange(1,numberSiO2+1):
        tempSiO2 = float(input('Enter the SiO2 thickness in nanometres: '))
        tempSiO2 = tempSiO2/1000
        if tempSiO2 < 0 or tempSiO2 > 500:
            print('Only wavelengths between 300nm and 500nm are valid, please try again')
        SiO2List.append(tempSiO2)
elif graphType in [3, 'three']:
    graphType = 3
    SiOrPDMS = raw_input("Are you looking at a flake on a Silicon or PDMS substrate? ")
    if SiOrPDMS.lower() in ['pdms','p','pd','pms','not si', 'nosi','two','second','other']:
        nvalsSi = nvalsPDMS
        nvalsSiO2 = nvalsPDMS
        setThickness = 300 #just so there is something here, really the PDMS is taken to be semi-infinite
    else:
        setThickness = float(input("Enter the thickness of the SiO2 substrate you would like to use in nanometres: "))
        setThickness = setThickness/1000
    numberThicknesses = input("Enter the number of {} thicknesses you would like to look at: ".format(formalName))
    thicknessList = []
    for i in xrange(1,numberThicknesses+1):
        tempThickness = float(input('Enter the {} thickness in nanometres: '.format(formalName)))
        tempThickness = tempThickness/1000
        if tempThickness < 0 or tempThickness > 300:
            print('Only thicknesses between 0 and 300 nanometres are valid, please try again')
        thicknessList.append(tempThickness)
elif graphType in [4, 'four']:
    graphType = 4
    setThickness = float(input("Enter the thickness of {}, in nanometres. A good idea is to choose a monolayer ({}nm): ".format(formalName,monoLayer*1000)))
    setThickness = setThickness/1000
elif graphType in [5, 'five']:
    graphType = 5
    SiOrPDMS = raw_input("Are you looking at a flake on a Silicon or PDMS substrate? ")
    if SiOrPDMS.lower() in ['pdms','p','pd','pms','not si', 'nosi','two','second','other']:
        nvalsSi = nvalsPDMS
        nvalsSiO2 = nvalsPDMS
        setThickness = 300 #just so there is something here, really the PDMS is taken to be semi-infinite
    else:
        setThickness = float(input("Enter the thickness of the SiO2, in nanometres: "))
    numberMonolayers = float(input("Enter the maximum number of monolayers you would like to examine: "))
    print("We will look at {} thicknesses from 0 to {} monolayers".format(formalName,numberMonolayers))
    setThickness = setThickness/1000
elif graphType in [6, 'six']:
    graphType = 6
    #print("We will look at the contrast for 0-50 monolayers of {}.".format(formalName))
    monolayThick = monoLayer
    
"""Plot the refractive index"""
if graphType == 1:
    print('Showing the refractive index plots of {}, Si, and SiO2'.format(formalName))
    if materialType == 'hbn':
        RefracIndexBN, axBN = plt.subplots(1,1)
        ldaBN = [i*1000 for i in ldaBN]
        l1 = axBN.plot(ldaBN,BNnreal, label = 'BN real')
        l2 = axBN.plot(ldaBN,BNnimag, label = 'BN imaginary')
        l3 = axBN.plot(lda*1000,nvalsBN_2, label = 'BN_2 real')
        plt.legend(); plt.xlim(400,800); plt.title('Boron Nitride refractive index plot')
        plt.xlabel('Wavelength (um)'); plt.ylabel('Refractive index')
    elif materialType == 'grp':
        print('\nThere are two refractive indices for graphene cited.')
        print('Blake et al. in Making graphene visible use n=2.6-1.3j')
        print('Another paper uses n=2.0-1.1j')
        print('The refractive index is currently set to n=2.6-1.3j however you can change this in the program\n')
    elif materialType == 'mica':
        print('Unknown refractive index, refractive index for a wavelength of 542nm has been estimated to be 1.19 using this program')
    elif materialType == 'sb':
        RefracIndexSb, axSb = plt.subplots(1,1)
        ldaSb = [i*1000 for i in ldaSb]
        l1 = axSb.plot(ldaSb,Sbnreal, label = 'Real component')
        l2 = axSb.plot(ldaSb,Sbnimag, label = 'Imaginary component')
        plt.legend(); plt.title('Antimony refractive index plot'); plt.xlim(0.4,0.8)
        plt.xlabel('Wavelength (um)'); plt.ylabel('Refractive index')
    elif materialType == 'nacl':
        RefracIndexNaCl, axNaCl = plt.subplots(1,1)
        l1 = axNaCl.plot(lda*1000,nvalsNaCl, label = 'Real component')
        plt.legend(); plt.xlim(400,800); plt.title('Sodium Chloride refractive index plot')
        plt.xlabel('Wavelength (um)'); plt.ylabel('Refractive index')
    elif materialType == 'wse2':
        RefracIndexWSe2, axWSe2 = plt.subplots(1,1)
        l1 = axWSe2.plot(lda*1000,nvalsWSe2.real, label = 'Real component')
        l2 = axWSe2.plot(lda*1000,-nvalsWSe2.imag, label = 'Imaginary component')
        plt.legend(); plt.xlim(400,800); plt.title('WSe2 refractive index plot')
        plt.xlabel('Wavelength (um)'); plt.ylabel('Refractive index')
    elif materialType == 'mos2':
        RefracIndexMoS2, axMoS2 = plt.subplots(1,1)
        l1 = axMoS2.plot(lda*1000,nvalsMoS2.real, label = 'Real component')
        l2 = axMoS2.plot(lda*1000,-nvalsMoS2.imag, label = 'Imaginary component')
        plt.legend(); plt.xlim(0.4,0.8); plt.title('MoS2 refractive index plot')
        plt.xlabel('Wavelength (um)'); plt.ylabel('Refractive index')
    elif materialType == 'nbse2':
        RefracIndexNbSe2, axNbSe2 = plt.subplots(1,1)
        l1 = axNbSe2.plot(lda*1000,nvalsNbSe2.real, label = 'Real component')
        l2 = axNbSe2.plot(lda*1000,-nvalsNbSe2.imag, label = 'Imaginary component')
        plt.legend(); plt.xlim(400,800); plt.title('NbSe2 refractive index plot')
        plt.xlabel('Wavelength (um)'); plt.ylabel('Refractive index')
    else:
        RefracIndex, ax = plt.subplots(1,1)
        l1 = ax.plot(lda*1000,nvals, label = 'Real component')
        plt.legend(); plt.xlim(400,800); plt.title('{} refractive index plot'.format(formalName))
        plt.xlabel('Wavelength (um)'); plt.ylabel('Refractive index')
    #Si and SiO2
    print('There are two good empirical sources for the refractive index of Si. Since they are very similar, only one is shown here, \
although both are available and you can choose which one you want to use in the program.')
    RefracIndexSi, axSi = plt.subplots(1,1)
    l1 = axSi.plot(lda,nvalsSi,label='Real component of Si')
    l2 = axSi.plot(lda,nvalsSi*(1j),label='Imaginary component of Si')
    l3 = axSi.plot(lda,nvalsSiO2,label='Silicon Oxide')
    plt.legend()
    plt.title('Refractive index of Silicon')
    plt.xlabel('Wavelength (micrometres)'); plt.ylabel('Refractive index')
   
"""Define the formula for calculating the contrast, from the paper"""
def rvals(n0,n1,n2,n3):
    r1 = (n0 - n1)/(n0 + n1)
    r2 = (n1 - n2)/(n1 + n2)
    r3 = (n2 - n3)/(n2 + n3)
    return r1, r2, r3
    
def phasevals(n1,n2,d1,d2,lda):
    phase1 = (2*(np.pi)*n1*d1)/lda
    phase2 = (2*(np.pi)*n2*d2)/lda
    return phase1,phase2
    
## write I as a function so we can insert value n1
def I(n0,n1,n2,n3,d1,d2,lda):
    r1, r2, r3 = rvals(n0,n1,n2,n3)
    phase1, phase2 = phasevals(n1,n2,d1,d2,lda)    
    Ipart1 = r1*np.exp(1j*(phase1+phase2)) + r2*np.exp(-1j*(phase1-phase2)) + \
             r3*np.exp(-1j*(phase1+phase2)) + r1*r2*r3*np.exp(1j*(phase1-phase2))
    Ipart2 = np.exp(1j*(phase1+phase2)) + r1*r2*np.exp(-1j*(phase1-phase2)) + \
             r1*r3*np.exp(-1j*(phase1+phase2)) + r2*r3*np.exp(1j*(phase1-phase2))
    Ipart2 = Ipart2**(-1)
    Itotal = abs(Ipart1*Ipart2)**2
    return Itotal

def Contrast(n0,n1,n2,n3, d1,d2,lda):
    #This contrast is known as the Weber contrast.
    C = (I(n0,1,n2,n3,d1,d2,lda)-I(n0,n1,n2,n3,d1,d2,lda))/I(n0,1,n2,n3,d1,d2,lda)
    """The Michelson contrast (below) is an alternative that is widely used for gratings.
    C = (I(n0,n1,n2,n3,d1,d2,lda)-I(n0,1,n2,n3,d1,d2,lda))/(I(n0,n1,n2,n3,d1,d2,lda)+I(n0,1,n2,n3,d1,d2,lda))"""
    return C

"""Set refractive indices from top material downwards"""
n0 = 1 #refractive index of air
#n1 is the refractive index of our material
n2 = nvalsSiO2 #set n2 = SiO2
n3 = nvalsSi #set n3 = Si

"""Reverse engineering to see which refractive index for mica works best"""
nvalsList = np.linspace(0.00,4.00,100)
def contrastForNVALS(nvalsList):
    contrastNVAL = np.zeros(len(nvalsList))
    for i in range(0,len(nvalsList)):
        #contrastNVAL[i] = Contrast(n0,nvalsList[i],nvalsSiO2[553],nvalsSi[553],0.035,0.3,0.610)
        contrastNVAL[i] = Contrast(n0,nvalsList[i],nvalsSiO2[710],nvalsSi[710],0.035,0.3,0.673) 
    return contrastNVAL

"""Calculate the contrast for different materials, as a function of SiO2 thickness"""    
def contrastForSiO2Thickness(numberSiO2, SiO2List, setThickness):
    contrastSiO2Thick = np.zeros((numberSiO2, len(lda)))
    for i in range(0, len(lda)):
        for j in range(0, numberSiO2):
            contrastSiO2Thick[j][i] = Contrast(n0,nvals[i],n2[i],n3[i],setThickness,SiO2List[j],lda[i]) 
    return contrastSiO2Thick
    
#plots the contrast
if graphType == 2:
    contrastSiO2Thick = contrastForSiO2Thickness(numberSiO2, SiO2List, setThickness)
    contrastSiO2plot, axSiO2Plot = plt.subplots(1,1)
    plt.title('Contrast plot for {}nm of {} for different SiO2 thickness'.format(setThickness*1000,formalName))
    for x in range(0, numberSiO2):
        globals()['l%s' % x] = axSiO2Plot.plot(lda*1000, contrastSiO2Thick[x], label = '{} nm'.format(SiO2List[x]*1000))
    plt.legend(); plt.xlabel('Wavelength (nanometres)'); plt.ylabel('Contrast')

"""Repeat this, but now varying the material thickness. First, calculate the contrast arrays"""
def contrastFor2DThickness(numberThicknesses, thicknessList, setThickness):
    contrast2DThick = np.zeros((numberThicknesses, len(lda)))
    for i in range(0, numberThicknesses):
        for j in range(0, len(lda)):
            contrast2DThick[i][j] = Contrast(n0,nvals[j],n2[j],n3[j],thicknessList[i],setThickness,lda[j]) 
    return contrast2DThick

"""Plot the contrast for a certain SiO2 thickness and for many material thicknesses"""
if graphType == 3:
    contrast2DThick = contrastFor2DThickness(numberThicknesses, thicknessList, setThickness)
    contrast2Dplot, ax2DPlot = plt.subplots(1,1)
    plt.title('Contrast plot for {} for different material thickness'.format(formalName))
    for x in range(0, numberThicknesses):
        globals()['l%s' % x] = ax2DPlot.plot(lda*1000, contrast2DThick[x], label = '{} nm'.format(thicknessList[x]*1000))
        globals()['inter%s' % (x+1)] = interpolate.interp1d(lda, contrast2DThick[x])
    plt.legend(); plt.xlabel('Wavelength (nanometres)'); plt.ylabel('Contrast')
    print("If you would like to know the contrast for a certain wavelength for a certain thickness, type: \n \
    interX(wavelength value) where X is the number of the thickness you would like to examine.")
    
#Create a plot of wavelength against SiO2 thickness with contrast as the colour
def ColourContrastForSiO2(SiO2Vals,lda):
    return Contrast(n0,nvals,SiO2realn(lda),(Sirealn(lda) + (Siimagn(lda)*(-1j))),setThickness,SiO2Vals,lda) 
#Create a plot of wavelength against monolayer thickness with contrast as the colour
def ColourContrastForMonolayers(matVals,lda):
    return Contrast(n0,nvals,SiO2realn(lda),(Sirealn(lda) + (Siimagn(lda)*(-1j))),matVals,setThickness,lda)

"""Create a colour contrast plot of wavelength against SiO2 thickness"""
if graphType == 4:
    SiO2Vals = np.linspace(0,350,1000) #0 to 0.3 in um
    lda = lda*1000
    for i in range(0,len(SiO2Vals)):
        SiO2Vals[i] = float(SiO2Vals[i])
    X,Y = meshgrid(SiO2Vals,lda) # grid of points
    X, Y = X[::-1],Y[::-1]
    Z = ColourContrastForSiO2(X/1000, Y/1000) # evaluation of the function on the grid
    im2 = plt.imshow(Z,cmap=cm.RdBu, vmin=-np.abs(Z).max(), vmax=np.abs(Z).max(), extent=[SiO2Vals.min(), SiO2Vals.max(), lda.min(), lda.max()], aspect='auto') # drawing the function
    plt.colorbar(im2) # adding the colobar on the right
    plt.clim(Z.min(), Z.max()) #adjust the range of the colorbar
    plt.title('Contrast plot for {}'.format(formalName)); plt.xlabel('SiO2 thickness'); plt.ylabel('Wavelength (nanometres)')    

"""Create a colour contrast plot of wavelength against material thickness in monolayers"""
if graphType == 5:
    matVals = np.linspace(0,numberMonolayers,1000) #0 to 0.3 in um
    lda = lda*1000
    X,Y = meshgrid(matVals,lda) # grid of points
    X, Y = X[::-1],Y[::-1]
    Z = ColourContrastForMonolayers((X*monoLayer), Y/1000) # evaluation of the function on the grid
    im3 = plt.imshow(Z,cmap=cm.RdBu, vmin=-np.abs(Z).max(), vmax=np.abs(Z).max(), extent=[matVals.min(), matVals.max(), lda.min(), lda.max()], aspect='auto') # drawing the function
    plt.colorbar(im3) # adding the colobar on the right
    plt.clim(Z.min(), Z.max())
    plt.title('Contrast plot for {}'.format(formalName)); plt.xlabel('{} monolayers'.format(formalName)); plt.ylabel('Wavelength (nanometres)')    

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

"""Create a plot for our image test,
530nm wavelength light
Create a contrast plot for many thicknesses
0.0034 is the smallest thickness (monolayer)"""
if graphType == 6:
    print('Select the image of your flake (may have to minimise this window)')
    execfile("ContrastMeasurement.py") #run the program 'Contrast from Image' within this program
    SiOrPDMS = raw_input("Are you looking at a flake on a Silicon or PDMS substrate? ")
    if SiOrPDMS.lower() in ['pdms','p','pd','pms','not si', 'nosi','two','second','other']:
        nvalsSi = nvalsPDMS
        nvalsSiO2 = nvalsPDMS
        SiO2setThick = 300./1000 #just so there is something here, really the PDMS is taken to be semi-infinite
    else:
        SiO2setThick = float(input("Enter the SiO2 thickness in nm: "))
        SiO2setThick = SiO2setThick/1000
    WavelsetThick = float(input("Enter the wavelength used in nm (we have filters for 542nm, 610n, 673nm): "))
    WavelsetThick = WavelsetThick/1000
    #now we need to convert that contrast to a thickness
    thickness = np.linspace(0,(monolayThick*100),101) #50 monolayers thickness
    NoThick = len(thickness)
    idx = find_nearest(lda,WavelsetThick) #find the index of the nearest wavelength in our lambda array to the wavelength specified by the user
    wlgth = WavelsetThick #shorten the name
    wlgthErr = 0.002 #let the error in wavelength due to bandpass filter be 2nm
    nValue = nvals[idx]
    #CForT denotes Contrast For Thickness
    CForT = np.zeros(NoThick) #SiO2 thickness entered
    CForTWErr1 = np.zeros(NoThick)  #SiO2 thickness entered plus Wavelength error
    CForTWErr2 = np.zeros(NoThick) #SiO2 thickness entered minus Wavelength error
    CForTP5 = np.zeros(NoThick) #SiO2 thickness with error
    CForTP5Err1 = np.zeros(NoThick) #SiO2 thickness with error plus wavelength error
    CForTP5Err2 = np.zeros(NoThick) #SiO2 thickness with error minus wavelength error
    CForTM5 = np.zeros(NoThick) #SiO2 thickness with error
    CForTM5Err1 = np.zeros(NoThick) #SiO2 thickness with error plus wavelength error
    CForTM5Err2 = np.zeros(NoThick) #SiO2 thickness with error minus wavelength error
    minArray = []
    maxArray = []
    if SiOrPDMS == 'pdms':
        for i in range(0,NoThick):
            PDMSnValue = nvalsPDMS[find_nearest(lda,wlgth)] #PDMS refractive index for the wavelength specified
            CForT[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick,wlgth)
            CForTWErr1[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick,wlgth+wlgthErr) #PDMS refractive index doesn't change with wavelength
            CForTWErr2[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick,wlgth-wlgthErr)
            CForTP5[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick*1.05,wlgth)
            CForTP5Err1[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick*1.05,wlgth+wlgthErr)
            CForTP5Err2[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick*1.05,wlgth-wlgthErr)
            CForTM5[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick*0.95,wlgth)
            CForTM5Err1[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick*0.95,wlgth+wlgthErr)
            CForTM5Err2[i] = Contrast(n0,nValue,PDMSnValue,PDMSnValue,thickness[i],SiO2setThick*0.95,wlgth-wlgthErr)
            minArray.append(np.min([ CForT[i],CForTP5[i],CForTM5[i],CForTWErr1[i],CForTWErr2[i],CForTP5Err1[i],CForTP5Err2[i],CForTM5Err1[i],CForTM5Err2[i] ]))
            maxArray.append(np.max([ CForT[i],CForTP5[i],CForTM5[i],CForTWErr1[i],CForTWErr1[i],CForTWErr2[i],CForTP5Err1[i],CForTP5Err2[i],CForTM5Err1[i],CForTM5Err2[i] ]))
    else:        
        for i in range(0,NoThick):
            CForT[i] = Contrast(n0,nValue,SiO2realn(wlgth),(Sirealn(wlgth)+(Siimagn(wlgth)*(-1j))),thickness[i],SiO2setThick,wlgth)
            CForTWErr1[i] = Contrast(n0,nValue,SiO2realn(wlgth+wlgthErr),(Sirealn(wlgth+wlgthErr)+(Siimagn(wlgth+wlgthErr)*(-1j))),thickness[i],SiO2setThick,wlgth+wlgthErr)
            CForTWErr2[i] = Contrast(n0,nValue,SiO2realn(wlgth-wlgthErr),(Sirealn(wlgth-wlgthErr)+(Siimagn(wlgth-wlgthErr)*(-1j))),thickness[i],SiO2setThick,wlgth-wlgthErr)
            CForTP5[i] = Contrast(n0,nValue,SiO2realn(wlgth),(Sirealn(wlgth)+(Siimagn(wlgth)*(-1j))),thickness[i],SiO2setThick*1.05,wlgth)
            CForTP5Err1[i] = Contrast(n0,nValue,SiO2realn(wlgth+wlgthErr),(Sirealn(wlgth+wlgthErr)+(Siimagn(wlgth+wlgthErr)*(-1j))),thickness[i],SiO2setThick*1.05,wlgth+wlgthErr)
            CForTP5Err2[i] = Contrast(n0,nValue,SiO2realn(wlgth-wlgthErr),(Sirealn(wlgth-wlgthErr)+(Siimagn(wlgth-wlgthErr)*(-1j))),thickness[i],SiO2setThick*1.05,wlgth-wlgthErr)
            CForTM5[i] = Contrast(n0,nValue,SiO2realn(wlgth),(Sirealn(wlgth)+(Siimagn(wlgth)*(-1j))),thickness[i],SiO2setThick*0.95,wlgth)
            CForTM5Err1[i] = Contrast(n0,nValue,SiO2realn(wlgth+wlgthErr),(Sirealn(wlgth+wlgthErr)+(Siimagn(wlgth+wlgthErr)*(-1j))),thickness[i],SiO2setThick*0.95,wlgth+wlgthErr)
            CForTM5Err2[i] = Contrast(n0,nValue,SiO2realn(wlgth-wlgthErr),(Sirealn(wlgth-wlgthErr)+(Siimagn(wlgth-wlgthErr)*(-1j))),thickness[i],SiO2setThick*0.95,wlgth-wlgthErr)
            minArray.append(np.min([ CForT[i],CForTP5[i],CForTM5[i],CForTWErr1[i],CForTWErr2[i],CForTP5Err1[i],CForTP5Err2[i],CForTM5Err1[i],CForTM5Err2[i] ]))
            maxArray.append(np.max([ CForT[i],CForTP5[i],CForTM5[i],CForTWErr1[i],CForTWErr1[i],CForTWErr2[i],CForTP5Err1[i],CForTP5Err2[i],CForTM5Err1[i],CForTM5Err2[i] ]))
    #This graph has repeated values on each side, which would mess with the contrast calculations
    #we will assume the graphene flake is on the thinner side
    """Max val index must be positive or negative!!!"""
    #only look at before the maximum
    CForT = np.array(CForT).tolist()
    CForTWErr1 = np.array(CForTWErr1).tolist()
    CForTWErr2 = np.array(CForTWErr2).tolist()
    """Find the turning points"""
    turningPoint = -1 #default turning point is the last index of the array, i.e. there is no turning point!
    thickInput = str(0) #also needs to be initialised
    for i in range(1, len(CForT)-1):
        if ((CForT[i-1]<CForT[i] and CForT[i+1]<CForT[i]) or (CForT[i-1]>CForT[i] and CForT[i+1]>CForT[i])): #if there is a turning point
            turningPoint = i #label that index
            print("There is a turning point in the contrast for {} at {}nm".format(formalName,thickness[i]*1000))
            thickInput = raw_input("Is your flake generally thinner or thicker than that?: ")
    ContrastThickness, axThickness = plt.subplots(1,1)
    l1 = axThickness.plot(range(0,NoThick),CForT,label='Contrast for 300nm SiO2, 542nm wavelength')
    l10 = axThickness.plot(range(0,NoThick),minArray, label='Maximum contrast for error on SiO2 and wavelength')
    l11 = axThickness.plot(range(0,NoThick),maxArray, label='Minimum contrast for error on SiO2 and wavelength')
    plt.title('Contrast against thickness for {}, wavelength: {}um, SiO2: {}nm'.format(formalName,WavelsetThick*1000,SiO2setThick*1000))
    axThickness.fill_between(range(0,NoThick), minArray, maxArray, facecolor='red', interpolate=True, alpha=0.2)
    plt.xlabel('Thickness (monolayers)'); plt.ylabel('Contrast');

    """Next task: save contrast array from another image to a file, which we can then load from here"""
    totalDir = os.path.join(scriptDir, "ContrastFiles/*")
    txtDir = os.path.join(scriptDir, "ContrastFiles/*.txt")
    allFiles = glob.glob(totalDir) #import all the contrast files
    txtFiles = glob.glob(txtDir) #import all the .txt contrast files (just contrast, not averaged)
    avgContrastFile = max(allFiles, key=os.path.getctime) #find the most recent avgContrast file
    contrastFile = max(txtFiles, key=os.path.getctime) #fine the most recent contrast file    
    
    avgContrast = [line.rstrip('\n') for line in open('{}'.format(avgContrastFile))] #strip the carriage return
    for i in range(0,len(avgContrast)):
        avgContrast[i] = float(avgContrast[i]) #convert to floats
    measContrast = [line.rstrip('\n') for line in open('{}'.format(contrastFile))]
    for i in range(0,len(measContrast)):
        measContrast[i] = float(measContrast[i]) #convert to floats
        
    precision = 1000 #how close you want to match contrast values with thickness values
    InterpContrast = np.zeros(precision)
    InterpContrastErr1 = np.zeros(precision)
    InterpContrastErr2 = np.zeros(precision)
    ContrToThickness = np.zeros(len(avgContrast))
    ContrToThicknessErr1 = np.zeros(len(avgContrast))
    ContrToThicknessErr2 = np.zeros(len(avgContrast))
    InterpNoAvgContrast = np.zeros(precision)
    NoAvgContrToThickness = np.zeros(len(measContrast))
    """Interpolate the contrast as a function of monolayer thickness plot, but only up until the turning point!"""
    if thickInput in ['thinner', 't','Thinner','THINNER','thin','Thin','THIN']: #if the flake is thinner, plot UP TO the turning point  
        InterpCFT = interpolate.interp1d(range(0,turningPoint+1), CForT[:turningPoint+1]) #interpolate the contrast
        InterpCFTErr1 = interpolate.interp1d(range(0,turningPoint+1), minArray[:turningPoint+1])
        InterpCFTErr2 = interpolate.interp1d(range(0,turningPoint+1), maxArray[:turningPoint+1])
        for x in range(0,precision): #find interpolated values for the contrast in the plotted graph of contrast vs thicknesses
            vThick= np.linspace(0,turningPoint,precision)
            InterpContrast[x] = InterpCFT(vThick[x])
            InterpContrastErr1[x] = InterpCFTErr1(vThick[x])   
            InterpContrastErr2[x] = InterpCFTErr2(vThick[x])   
        for i in range(0,len(avgContrast)):
            idx = find_nearest(InterpContrast,avgContrast[i]) #find the contrast value closest to that in the array
            idxErr1 = find_nearest(InterpContrastErr1,avgContrast[i]) #find the contrast value closest to that in the array
            idxErr2 = find_nearest(InterpContrastErr2,avgContrast[i]) #find the contrast value closest to that in the array
            ContrToThickness[i] = vThick[idx] #then we need to match that to a thickness in the InterpContrast
            ContrToThicknessErr1[i] = vThick[idxErr1] #then we need to match that to a thickness in the InterpContrast
            ContrToThicknessErr2[i] = vThick[idxErr2] #then we need to match that to a thickness in the InterpContrast
        for i in range(0,len(measContrast)):
            noAvgIdx = find_nearest(InterpContrast,measContrast[i])
            NoAvgContrToThickness[i] = vThick[noAvgIdx]
    """If the user says the flake is on the thicker side, interpolate after the turning point"""
    if thickInput in ['thicker', 'T','Thicker','THICKER','thick','Thick','THICK','thck','thk','tick']: #if the flake is thick, plot AFTER the turning point  
        InterpCFT = interpolate.interp1d(range(turningPoint,len(CForT)), CForT[turningPoint:]) #interpolate the contrast
        InterpCFTErr1 = interpolate.interp1d(range(turningPoint,len(CForT)), minArray[turningPoint:])
        InterpCFTErr2 = interpolate.interp1d(range(turningPoint,len(CForT)), maxArray[turningPoint:])  
        for x in range(0,precision): #find interpolated values for the contrast in the plotted graph of contrast vs thicknesses
            vThick= np.linspace(turningPoint,len(CForT)-1,precision)
            InterpContrast[x] = InterpCFT(vThick[x])
            InterpContrastErr1[x] = InterpCFTErr1(vThick[x])   
            InterpContrastErr2[x] = InterpCFTErr2(vThick[x])  
        for i in range(0,len(avgContrast)):
            idx = find_nearest(InterpContrast,avgContrast[i]) #find the contrast value closest to that in the array
            idxErr1 = find_nearest(InterpContrastErr1,avgContrast[i]) #find the contrast value closest to that in the array
            idxErr2 = find_nearest(InterpContrastErr2,avgContrast[i]) #find the contrast value closest to that in the array
            ContrToThickness[i] = vThick[idx] #then we need to match that to a thickness in the InterpContrast
            ContrToThicknessErr1[i] = vThick[idxErr1] #then we need to match that to a thickness in the InterpContrast
            ContrToThicknessErr2[i] = vThick[idxErr2] #then we need to match that to a thickness in the InterpContrast
        for i in range(0,len(measContrast)):
            noAvgIdx = find_nearest(InterpContrast,measContrast[i])
            NoAvgContrToThickness[i] = vThick[noAvgIdx]
    else: pass
    """In the case that there is no turning point"""        
    if turningPoint == -1 or thickInput == str(0): #i.e. turning point is the last point on the array             
        InterpCFT = interpolate.interp1d(range(0,len(CForT)), CForT) #interpolate the contrast
        InterpCFTErr1 = interpolate.interp1d(range(0,len(CForT)), minArray)
        InterpCFTErr2 = interpolate.interp1d(range(0,len(CForT)), maxArray)  
        for x in range(0,precision): #find interpolated values for the contrast in the plotted graph of contrast vs thicknesses
            vThick= np.linspace(0,len(CForT)-1,precision)
            InterpContrast[x] = InterpCFT(vThick[x])
            InterpContrastErr1[x] = InterpCFTErr1(vThick[x])   
            InterpContrastErr2[x] = InterpCFTErr2(vThick[x])  
        for i in range(0,len(avgContrast)):
            idx = find_nearest(InterpContrast,avgContrast[i]) #find the contrast value closest to that in the array
            idxErr1 = find_nearest(InterpContrastErr1,avgContrast[i]) #find the contrast value closest to that in the array
            idxErr2 = find_nearest(InterpContrastErr2,avgContrast[i]) #find the contrast value closest to that in the array
            ContrToThickness[i] = vThick[idx] #then we need to match that to a thickness in the InterpContrast
            ContrToThicknessErr1[i] = vThick[idxErr1] #then we need to match that to a thickness in the InterpContrast
            ContrToThicknessErr2[i] = vThick[idxErr2] #then we need to match that to a thickness in the InterpContrast    
        for i in range(0,len(measContrast)):
            noAvgIdx = find_nearest(InterpContrast,measContrast[i])
            NoAvgContrToThickness[i] = vThick[noAvgIdx]
    """Plot the thickness for the line drawn in the image"""
    ThicknessPlot, axThick = plt.subplots(1,1)
    averaging = len(NoAvgContrToThickness)/len(ContrToThickness)
    l1 = axThick.plot(np.linspace(0,len(NoAvgContrToThickness),len(ContrToThickness)),ContrToThickness, label = '%s point average' %(avgNo))
    #l2 = axThick.plot(range(0,len(ContrToThickness)),ContrToThicknessErr1, alpha = 0.2)
    #l3 = axThick.plot(range(0,len(ContrToThickness)),ContrToThicknessErr2, alpha = 0.2)
    l4 = axThick.plot(range(0,len(NoAvgContrToThickness)),NoAvgContrToThickness,label = 'No average')
    plt.title('Thickness for {}, wavelength: {}um, SiO2: {}nm'.format(formalName,WavelsetThick*1000,SiO2setThick*1000))
    #axThick.fill_between(range(0,len(ContrToThickness)), ContrToThicknessErr2, ContrToThicknessErr1, facecolor='red', interpolate=True, alpha=0.2)
    plt.xlabel('Displacement along line'); plt.ylabel('Thickness (monolayers)'); plt.legend()

