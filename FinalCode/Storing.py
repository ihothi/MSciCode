import numpy as np
from astropy.io import fits
import os
from Project_Functions import MLAData,classification, Object,storing,StandardRebin

# Pixel Rejection Threshold
reject_no = 500
# Bin Size
Bin_Size=10
slash =  os.path.normpath("/")
# Plate data path
Platedir = os.path.normpath(slash+"share"+slash+"data1"+slash+"boss_data"+slash+"sas"+slash+"dr12"+slash+"boss"+slash+"spectro"+slash+"redux"+slash+"v5_7_2")
# What is the path for the data to be stored?
Bin_platedir =  os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"New_boss_data_rebinned_x30")
#List of all plates
plate_name = os.listdir(Platedir)
file="Full_Name.txt"
p = open(file, 'w')
#creating a file of all plates
for i in plate_name:
    p.write(i +'\n')

p.close()
#Putting plate names into an array
with open(file) as f:
    Files = f.read().splitlines()

print(np.str(len(Files)))
Spectra_Files=[]
i=0

while i< len(Files):
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
#Opening Files and Storing Data
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
#Opening Superset Data   
list = fits.open('Superset_DR12Q.fits',memmap=True)#opening file
supers=list[1].data # storing  BINTABLE extension data

#Storing the Superset Data
Full_Data = storing(PLATEIDs,supers)
# Matching the plate data to the superset to get classification
# As well as putting it into a suitable form for the 
#Machine learning algorithm (MLA)
X,Y,Train_z, Train_mag,And, In, wavst, ID, MJ, MatchedPlates = MLAData(Full_Data,BinInfos,Flux, log_wavst,ANDMASK,INVAR,MJDs,PLATEIDs)
plate_no = 0
all_Testwav=[]
wav_ratio = 10**0.0001


#Creating the Standard Wavelength Grid
while plate_no < len(wavst): 
    append_count=0
    cent_wav = 10**wavst[plate_no]
    wavelengths = []
    wavelengths.append(cent_wav)
    current_wav = cent_wav
    while append_count < (4600-1):
        current_wav = current_wav*wav_ratio
        wavelengths.append(current_wav)
        append_count=append_count+1
    plate_no=plate_no+1
    all_Testwav.append(wavelengths)
    
plate_no = 0
Xbin_plate=[]
true_Y =[]
true_Train_z=[]
true_ID =[]
true_MJ = []
true_MatchedPlates=[]
w_plate=[]
l_plate=[]
minimum = []
#Binning the pixel-flux for each object
while plate_no<len(X):
    plateX = X[plate_no]
    plateY = Y[plate_no]
    z0 = Train_z[plate_no]
    plateMask = And[plate_no]
    plateInv = In[plate_no]
    wave = all_Testwav[plate_no]
    pid  = ID[plate_no]  
    r,w, r_lambda,m,rebin_class,rebin_z,rebin_id = StandardRebin(plateX,wave,plateMask,plateInv ,Bin_Size,reject_no,plateY,z0,pid)
    true_Y.append(rebin_class)
    true_Train_z.append(rebin_z)
    true_ID.append(rebin_id)
    minimum.append(m)
    Xbin_plate.append(r)
    w_plate.append(w)
    l_plate.append(r_lambda)
    plate_no=plate_no+1
minbin = 500000
#Finding the Minimum number of bins 
for f in minimum:
    if f< minbin:
        minbin=f

i=0


#Storing the Plate files
First_Run=True
while i <len(MatchedPlates):
    C_Plate = MatchedPlates[i]
    a1 = Xbin_plate[i] 
    a2 = np.array(true_Y[i])
    a3=l_plate[i]
    a4=w_plate[i]
    a5 = true_Train_z[i]
    a6= true_ID[i]
    a8 = MJ[i]
    a7 = wavst[i]
    col1 = fits.Column(name='Bin_Flux', format='PD()', array=np.array(a1,dtype=np.object))
    col2 = fits.Column(name='Class', format='I', array=np.array(a2))
    col3 = fits.Column(name='Bin_Cent', format='PD()', array=np.array(a3,dtype=np.object))
    col4 = fits.Column(name='INVAR', format='PD()', array=np.array(a4,dtype=np.object))
    col5 = fits.Column(name='Redshift', format='D', array=np.array(a5))
    col6 = fits.Column(name='Name', format='20A', array=np.array(a6))
    cols = fits.ColDefs([col1, col2,col3, col4,col5,col6])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr['Plate'] = C_Plate
    prihdr['LogWav'] = a7
    prihdr['MJDs'] = a8
    prihdr['Min'] = minbin
    prihdu = fits.PrimaryHDU(header=prihdr)
    file_dir = Bin_platedir+slash+np.str(C_Plate)+slash
    try:   
        os.mkdir(file_dir)
        file_name = file_dir+np.str(C_Plate)+"_"+np.str(Bin_Size)+"_"+np.str(a8)+".fits"
        print(file_name)
        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(file_name)
        i=i+1
    # As a plate can be observed multiple times, store in the same folder
    # Just different MJD
    except FileExistsError as F:
        file_name = file_dir+np.str(C_Plate)+"_"+np.str(Bin_Size)+"_"+np.str(a8)+".fits"
        print(file_name)
        thdulist = fits.HDUList([prihdu, tbhdu])
        thdulist.writeto(file_name)
        i=i+1
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    