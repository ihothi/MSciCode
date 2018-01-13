import numpy as np
from astropy.io import fits
import os
from ProjectF import MLAData,classification, Object,storing







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
            primhdu_ = plate_[0]
            PLATEIDs.append(primhdu_.header['PLATEID'])
            ANDMASK.append( plate_[2].data)

        
list = fits.open('Superset_DR12Q.fits',memmap=True)#opening file
supers=list[1].data # storing  BINTABLE extension data
print("Restoring Data")
Full_Data = storing(PLATEIDs,supers)
X,Y,Train_z, Train_mag,And, In, wavst, ID, MJ, MatchedPlates = MLAData(Full_Data,BinInfos,Flux, log_wavst,ANDMASK,INVAR,MJDs,PLATEIDs)



    

i=0

i=0
print("Saving Files")
RejNoPlate = [] 
while i <len(MatchedPlates):
    no=0
    a7 = And[i]
    RejNo =0
    while no <len(a7):
        if a7[no] >0:
            RejNo = RejNo +1
        no=no+1
    RejNoPlate.append(RejNo) 
    i=i+1
        
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    