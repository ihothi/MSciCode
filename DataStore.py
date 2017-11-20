import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
from ProjectF import MLADataTest,classification, Object,storing
import random

## Loading in data
PlateDir = os.path.normpath("D:\Data\Plate_Name_Reduced.txt")
with open(PlateDir) as f:
    Spectra_Files = f.read().splitlines()
print('Yes')  
PLATEIDs = []
BinInfos = []
Flux = []
MJDs = []
log_wavst=[]
TrainingDir = os.path.normpath("D:\Data")
slash =  os.path.normpath("\\")
TrainingFolder =  os.path.normpath("\Training")
slash =  os.path.normpath("\\")
for spectrum in Spectra_Files:
    plate_ = fits.open( TrainingDir +TrainingFolder+slash+ spectrum ,memmap=True)
    Bin_info_ = plate_[5].data
    Flux_ = plate_[0].data
    primhdu_ = plate_[0]
    PLATEIDs.append(primhdu_.header['PLATEID'])
    log_wavst.append(primhdu_.header['COEFF0'])
    MJDs.append(primhdu_.header['MJD'])
    BinInfos.append(Bin_info_)
    Flux.append(Flux_)
    
list = fits.open(TrainingDir+slash+'Superset_DR12Q.fits',memmap=True)#opening file
print(len(log_wavst))
supers=list[1].data # storing  BINTABLE extension data

Full_Data = storing(PLATEIDs,supers)


X,Y,Train_z, Train_mag,wavst = MLADataTest(Full_Data,BinInfos,Flux, log_wavst)

i=0
while i <len(PLATEIDs):
    C_Plate = PLATEIDs[i]
    a1 = X[i] 
    a2 = Y[i]

    col1 = fits.Column(name='Bin_Flux', format='E', array=a1)
    col2 = fits.Column(name='Class', format='E', array=a2)
    cols = fits.ColDefs([col1, col2])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(TrainingDir+slash+"rebinned"+C_Plate+'.fits')
    i=i+1