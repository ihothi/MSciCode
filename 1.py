import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
from ProjectF import MLAData,classification, Object,storing,MLADataBin
import random
slash =  os.path.normpath("\\")
Platedir = os.path.normpath("D:"+slash+"share"+slash+"data1"+slash+"boss_data"+slash+"sas"+slash+"dr12"+slash+"boss"+slash+"spectro"+slash+"redux"+slash+"v5_7_0")

plate_name = os.listdir(Platedir)
file="FullPlate_Name.txt"
p = open(file, 'w')
for i in plate_name:
    print(i)
    p.write(i +'\n')

p.close()

with open(file) as f:
    Spectra_Files = f.read().splitlines() 
PLATEIDs = []
BinInfos = []
Flux = []
MJDs = []
log_wavst=[]
ORMASK=[]
ANDMASK=[]
INVAR=[]
print("Opening Files")
for f in Spectra_Files:
    file_list = os.listdir(Platedir+slash+f)
    for l in file_list:
        if 'spPlate' in l and ".fits"in l: 
            c=Platedir+slash+f+slash+l
            print(c)
            plate_ = fits.open(c,memmap=True)
            Bin_info_ = plate_[5].data
            Flux_ = plate_[0].data
            primhdu_ = plate_[0]
            PLATEIDs.append(primhdu_.header['PLATEID'])
            ORMASK.append( plate_[3].data)
            ANDMASK.append( plate_[2].data)
            INVAR.append( plate_[1].data)
            log_wavst.append(primhdu_.header['COEFF0'])
            MJDs.append(primhdu_.header['MJD'])
            BinInfos.append(Bin_info_)
            Flux.append(Flux_)
        
list = fits.open('Superset_DR12Q.fits',memmap=True)#opening file
supers=list[1].data # storing  BINTABLE extension data
print("Restoring Data")
Full_Data = storing(PLATEIDs,supers)
X,Y,Train_z, Train_mag,And, In, wavst, ID = MLAData(Full_Data,BinInfos,Flux, log_wavst,ANDMASK,INVAR)

i=0
print("Saving Files")
while i <len(PLATEIDs):
    C_Plate = PLATEIDs[i]
    a1 = X[i] 
    a2 = np.array(Y[i])
    if len(a1)==0 & len(a2)==0:
        i=i+1
    else: 
        a3=And[i]
        a4=In[i]
        a5 = Train_z[i]
        a6=ID[i]
        a7 = wavst[i]
        col1 = fits.Column(name='Bin_Flux', format='PD()', array=np.array(a1,dtype=np.object))
        col2 = fits.Column(name='Class', format='I', array=np.array(a2))
        col3 = fits.Column(name='ANDMASK', format='PD()', array=np.array(a3,dtype=np.object))
        col4 = fits.Column(name='INVAR', format='PD()', array=np.array(a4,dtype=np.object))
        col5 = fits.Column(name='Redshift', format='D', array=np.array(a5))
        col6 = fits.Column(name='Name', format='20A', array=np.array(a6))
        cols = fits.ColDefs([col1, col2,col3, col4,col5,col6])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        prihdr = fits.Header()
        prihdr['Plate'] = C_Plate
        prihdr['LogWav'] = a7
        prihdu = fits.PrimaryHDU(header=prihdr)
        file_name = "restore"+"/"+np.str(C_Plate)+'.fits'
        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(file_name)
        i=i+1