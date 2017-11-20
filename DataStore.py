import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
from ProjectF import MLAData,classification, Object,storing
import random

## Loading in data
PlateDir = os.path.normpath("D:\Data\Plate_Name_Reduced.txt")
with open(PlateDir) as f:
    Spectra_Files = f.read().splitlines()
    
PLATEIDs = []
BinInfos = []
Flux = []
MJDs = []
log_wavst=[]
TrainingDir = os.path.normpath("D:\Data")
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

supers=list[1].data # storing  BINTABLE extension data
    