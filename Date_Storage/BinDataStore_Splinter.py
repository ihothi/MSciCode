import numpy as np
from astropy.io import fits
import os
from ProjectF import MLAData,classification, Object,storing

def StandardRebin(plateX, wavelength,ANDMASK,INVAR ,Bin_Size ):
    
    rebin = []
    rebin_weight=[]
    minimum=100000000000
    bin_wav = Bin_Size*0.829026074968624
    wavelength=wavelength[:4600]
    rebinwav=[]
 
    obj=0
    while obj < len(plateX):
        current_flux = plateX[obj]
        current_mask = ANDMASK[obj]
        current_inv = INVAR[obj]
        rebin_flux=[]
        weight_=[]
        cwav=[]
        pix=0
        First = True
        bin_max = 3600
        while pix < len(current_flux):
            bin_no=0
            x_flux=0
            W=0
            wav=0
            if First: 
                cwav.append(bin_max+(bin_wav*0.5))
                bin_max = bin_max+bin_wav
                
                First = False
            else: 
                bin_max = bin_max+bin_wav
                cwav.append(bin_max+(bin_wav*0.5))
            while wav<bin_max:
                
                if pix< len(current_flux):
                    w_ = 0
                    m = current_mask[pix]
                    wav = wavelength[pix]
                    if m>0:
                        w_=0
                    else:
                        w_ = current_inv[pix]
                    x_flux = x_flux +(current_flux[pix]*w_)
                    W=W+w_
                    pix=pix+1
                else: 
                    pix=pix+1
                    wav=bin_max
                    
    
            if W == 0:
                x_flux=0
            else: 
                x_flux = x_flux/W
            rebin_flux.append(x_flux)

            weight_.append(W)
        rebinwav.append(cwav)
        rebin.append(rebin_flux)
        rebin_weight.append(weight_)
        if len(rebin_flux)<minimum:
            minimum=len(rebin_flux)
        obj=obj+1
    return rebin, rebin_weight,rebinwav,minimum





Bin_Size=10#np.int(input("Please Enter bin size: "))
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

print(np.str(len(Files)))
Spectra_Files=[]
i=1601
while i< 2442:
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
plate_no = 0
all_Testwav=[]
wav_ratio = 10**0.0001

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
w_plate=[]
l_plate=[]
minimum =0
print('Re-binning Test')
while plate_no<len(X):
    plateX = X[plate_no]
    plateMask = And[plate_no]
    plateInv = In[plate_no]
    wave = all_Testwav[plate_no]
    r,w, r_lambda,minimum = StandardRebin(plateX,wave,plateMask,plateInv ,Bin_Size )
    Xbin_plate.append(r)
    w_plate.append(w)
    l_plate.append(r_lambda)
    plate_no=plate_no+1
        
i=0

i=0
print("Saving Files")
First_Run=True
while i <len(MatchedPlates):
    C_Plate = MatchedPlates[i]
    a1 = Xbin_plate[i] 
    a2 = np.array(Y[i])
    if len(a1)==0 & len(a2)==0:
        i=i+1
    else: 
        a3=l_plate[i]
        a4=w_plate[i]
        a5 = Train_z[i]
        a6= ID[i]
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
        prihdr['Min'] = minimum
        prihdu = fits.PrimaryHDU(header=prihdr)
        file_dir = Bin_platedir+slash+np.str(C_Plate)+slash
        try:   
            os.mkdir(file_dir)
            file_name = file_dir+np.str(C_Plate)+"_"+np.str(Bin_Size)+"_"+np.str(a8)+".fits"
            print(file_name)
            thdulist = fits.HDUList([prihdu, tbhdu])
            thdulist.writeto(file_name)
            i=i+1
        except FileExistsError as F:
            file_name = file_dir+np.str(C_Plate)+"_"+np.str(Bin_Size)+"_"+np.str(a8)+".fits"
            print(file_name)
            thdulist = fits.HDUList([prihdu, tbhdu])
            thdulist.writeto(file_name)
            i=i+1
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    