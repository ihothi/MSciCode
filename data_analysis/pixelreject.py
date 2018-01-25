import numpy as np
from astropy.io import fits
import os
from ProjectF import MLAData,classification, Object,storing
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt







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
i=0
while i< 500:
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



    

i=0

i=0
print("Counting")
RejNoPlate = [] 
while i <len(MatchedPlates):
    no=0
    a7 = And[i]
    RejNo =0
    while no <len(a7):
        obj_and = a7[no]
        pixel_count = 0
        RejNo =0
        while pixel_count < len(obj_and):
            if obj_and[pixel_count] >0:
                RejNo = RejNo +1
            pixel_count=pixel_count+1
        no=no+1
        RejNoPlate.append(RejNo) 
    i=i+1

print(len(RejNoPlate))
print(RejNoPlate[0])
       
    
print("plotting")    
plt.hist(RejNoPlate)
plt.title('Rejection Plots')
plt.xlabel('Number of rejected pixels')
plt.ylabel('No. of Objects')
plt.savefig("Pixel_Rejection.png")
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    