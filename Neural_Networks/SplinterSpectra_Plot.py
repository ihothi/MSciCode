import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
from ProjectF import MLAData,classification, Object,storing,MLADataBin
from random import randint






Bin_Size = 10#np.int(input("Please Enter bin size: "))
slash =  os.path.normpath("/")
Platedir = os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"boss_data_rebinned_x10"+slash)
plate_name = os.listdir(Platedir)
file="FullPlate_Name.txt"
p = open(file, 'w')
for i in plate_name:
    p.write(i +'\n')

p.close()


with open(file) as f:
    FilesBin = f.read().splitlines() 
    
No_TrainPlates = 50 ##Note there are only 1200 plates total      
i=0
Spectra_BinFiles=[]
while i<No_TrainPlates :
    Spectra_BinFiles.append(FilesBin[i])
    i=i+1

wav_logBin=[]
Full_tableBin = []
PLATEIDsBin = []
print("Opening Files")
for f in Spectra_BinFiles:
    file_list = os.listdir(Platedir+slash+f)
    for l in file_list:
        c=Platedir+slash+f+slash+l
        plate_ = fits.open(c,memmap=True)
        primhdu_ = plate_[0]
        hdu =  plate_[1].data
        Full_tableBin.append(hdu)
        PLATEIDsBin.append(primhdu_.header['Plate'])
        wav_logBin.append(primhdu_.header['LogWav'])
        
plate_r  = []
plate_n  = []
plate_c  = []
plate_x  = []
plate_y  = []
plate_id =[]
wavelengths =[]
i = 0
while i < len(PLATEIDsBin):
    plate_hdu  = Full_tableBin[i]
    obj = 0
    while obj <len(plate_hdu):
        currentobj = plate_hdu[obj]
        Currentplate_z = currentobj[4]
        Currentplate_y = currentobj[1]
        Currentplate_x = currentobj[0]
        Currentplate_n = currentobj[5]
        wavelengths_plate = currentobj[2]
        plate_r.append(Currentplate_z)
        plate_n.append(Currentplate_n)
        plate_x.append(Currentplate_x[:800])
        plate_y.append(Currentplate_y)
        wavelengths.append(wavelengths_plate)
        plate_id.append(PLATEIDsBin[i])
        obj = obj +1
    i = i+1
    
    
    
    
slash =  os.path.normpath("/")
Platedir = os.path.normpath(slash+"share"+slash+"data1"+slash+"boss_data"+slash+"sas"+slash+"dr12"+slash+"boss"+slash+"spectro"+slash+"redux"+slash+"v5_7_0")
Bin_platedir =  os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"boss_data_rebinned_x10")
plate_name = os.listdir(Platedir)
file="Full_Name.txt"
p = open(file, 'w')
for i in plate_name:
    p.write(i +'\n')

p.close()

with open(file) as f:
    Files = f.read().splitlines()


Spectra_Files=[]
i=0
while i< 250:
    Spectra_Files.append(Files[i])
    i=i+1
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
X,Y,Train_z, Train_mag,And, In, wavst, ID, MJ, MatchedPlates = MLAData(Full_Data,BinInfos,Flux, log_wavst,ANDMASK,INVAR,MJDs,PLATEIDs)

X = X[:50]
Y= Y[:50]

print(np.str(len(X)))

Y_Full=[]
X_Full=[]  
plates_no =0
WFull = []
while plates_no < len(X):
    X_plate = X[plates_no] 
    Y_plate =Y[plates_no]
    plate_w=wavst[plates_no]
    obj_no = 0
    while obj_no < len(X_plate):
        X_Full.append(X_plate[obj_no])
        Y_Full.append(Y_plate[obj_no])
        WFull.append(plate_w)
        obj_no=obj_no+1
    plates_no= plates_no+1
    
  

not_found = True

X_Bin = []
X_plot = []
Bin_Waves = []
Plot_st = 0

no=0
if not_found:
     no = randint(0,len(plate_y))
     class_ = plate_y[no]
     
     if class_ == 1 or class_ == 4:
         print("Class is: ",np.str(class_))
         not_found= False


print(len(plate_x))
print(len(X_Full))
X_Bin = plate_x[no]
X_plot = X_Full[no]
print(np.str(X_plot))
Plot_st = WFull[no]
        
        

i=0
wav_ratio = 10**0.0001 
bin_wav = (10*0.829026074968624)
binmax = 3600*bin_wav
binmin = 3600

while i<800:
    Bin_Waves.append(binmin+(bin_wav*0.5))
    binmin = binmax
    binmax = binmin+bin_wav
    i=i+1
    
              
        
wave_ = []
total_pixel = len(X_plot)
cent_wav = 10**Plot_st
wave_.append(cent_wav)

pixel_count = 0



pixel_count=0
while pixel_count<(total_pixel-1):
    next_wave = cent_wav*wav_ratio
    wave_.append(next_wave)
    pixel_count=pixel_count+1


print("Plotting")
print(len(wave_))
print(len(X_plot))
plt.plot(wave_,X_plot)
plt.plot(Bin_Waves,X_Bin)
plt.ylabel('Flux [$10^{-17}$ $erg/s/cm^2/$'r'$\AA$'']')
plt.xlabel('Wavelength' r' $\AA$')
plt.title('Quasar spectra overplot')
plt.savefig('quasarbinplot.png')
    
    
        
        
    
    
    
    
    
