import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import os
from ProjectF import MLAData,classification, Object,storing,Rebin,StandardRebin
import random

slash =  os.path.normpath("/")
Platedir = os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"boss_data_rebinned_x10"+slash)
plate_name = os.listdir(Platedir)
file="FullPlate_Name.txt"
p = open(file, 'w')
for i in plate_name:
    p.write(i +'\n')

p.close()

with open(file) as f:
    Files = f.read().splitlines() 
    
i=0
print(len(Files))
Spectra_Files=[]
while i< 1200:
    if i == 64:
        print(Files[i])
        i=i+1
        
        
    else:
        Spectra_Files.append(Files[i])
        i=i+1


wav_log=[]
Full_table = []
PLATEIDs = []
print("Opening Files")
for f in Spectra_Files:
    file_list = os.listdir(Platedir+slash+f)
    for l in file_list:
        c=Platedir+slash+f+slash+l
        print(c)
        plate_ = fits.open(c,memmap=True)
        primhdu_ = plate_[0]
        hdu =  plate_[1].data
        Full_table.append(hdu)
        PLATEIDs.append(primhdu_.header['Plate'])
        wav_log.append(primhdu_.header['LogWav'])
        

plate_and = []
i = 0
while i < len(PLATEIDs):
    plate_hdu  = Full_table[i]
    obj = 0
    while obj <len(plate_hdu):
        currentobj = plate_hdu[obj]
        Currentplate_and = currentobj[4]
        plate_and.append(Currentplate_and)
        obj = obj +1
    i = i+1
    
wav_ratio = 10**0.0001
no = 0
Rej=[]

while no< len(plate_and):
    rej_count =0
    obj = plate_and[no]
    pixel =0
    while pixel < len(obj):
        current_pixel = obj[pixel]
        if current_pixel > 0:
            rej_count = rej_count+1
        pixel = pixel + 1
    no=no+1
   
plt.hist(Rej)
plt.title('Rejection Plots')
plt.xlabel('Number of rejected pixels')
plt.ylabel('No. of Objects')
plt.savefig("Pixel Rejection.png")