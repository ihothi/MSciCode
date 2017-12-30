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
from random import randint

def StandardRebin(plateX, wavelength,ANDMASK,INVAR ,Bin_Size ):
    
    rebin = []
    rebin_weight=[]
    bin_wav = Bin_Size*0.829026074968624
    wavelength=wavelength[:4600]
    rebinwav=[]
 
    obj=0
    while obj < len(plateX):
        last_wav=0
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
        obj=obj+1
        print("Rebin Progress: "+ np.str((obj*100)/len(plateX)))
    return rebin, rebin_weight,rebinwav














Bin_Size = 10#np.int(input("Please Enter bin size: "))
slash =  os.path.normpath("/")
Platedir = os.path.normpath(slash+"share"+slash+"data1"+slash+"boss_data"+slash+"sas"+slash+"dr12"+slash+"boss"+slash+"spectro"+slash+"redux"+slash+"v5_7_0")
plate_name = os.listdir(Platedir)
file="FullPlate_Name.txt"
p = open(file, 'w')
for i in plate_name:
    p.write(i +'\n')

p.close()

with open(file) as f:
    Files = f.read().splitlines() 
    
i=0
Spectra_Files=[]
while i< 500:
    a = randint(0, 2300)
    Spectra_Files.append(Files[a])
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
wav_ratio = 10**0.0001
plate_no = 0
all_wav=[]

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
    all_wav.append(wavelengths)
    
plate_no = 0
X_plate=[]
print('Re-binning')
while plate_no<len(X):
    plateX = X[plate_no]
    plateMask = And[plate_no]
    plateInv = In[plate_no]
    wave = all_wav[plate_no]
    r,w, r_weight = StandardRebin(plateX,wave,plateMask,plateInv ,Bin_Size )
    X_plate.append(r)
    plate_no=plate_no+1
    
X_Full = []
Y_Full = []
p = 0

while p < len(X_plate):
    CurrentplateX = X_plate[p]
    CurrentplateY=Y[p]
    n=0
    while n<len(CurrentplateX):
        X_Full.append(CurrentplateX[n][:800])
        Y_Full.append(CurrentplateY[n])
        n=n+1
    p=p+1
#Train_No = len(X_Full)/4
    
    





X_Test=[]
Y_Test = []
#i=0
#while i< Train_No:
#    a = randint(1, 100)
#    X_Test.append(X_Full[a])
#    del X_Full[a]
#    Y_Test.append(Y_Full[a])
#    del Y_Full[a]
#    i=i+1

i=0
Test_Files=[]
while i< 50:
    a = randint(0, 2300)
    Test_Files.append(Files[a])
    i=i+1
TestPLATEIDs = []
TestBinInfos = []
TestFlux = []
TestMJDs = []
Testlog_wavst=[]
TestORMASK=[]
TestANDMASK=[]
TestINVAR=[]
print("Opening Files")
for f in Test_Files:
    file_list = os.listdir(Platedir+slash+f)
    for l in file_list:
        if 'spPlate' in l and ".fits"in l: 
            c=Platedir+slash+f+slash+l
            print(c)
            plate_ = fits.open(c,memmap=True)
            Bin_info_ = plate_[5].data
            Flux_ = plate_[0].data
            primhdu_ = plate_[0]
            TestPLATEIDs.append(primhdu_.header['PLATEID'])
            TestORMASK.append( plate_[3].data)
            TestANDMASK.append( plate_[2].data)
            TestINVAR.append( plate_[1].data)
            Testlog_wavst.append(primhdu_.header['COEFF0'])
            TestMJDs.append(primhdu_.header['MJD'])
            TestBinInfos.append(Bin_info_)
            TestFlux.append(Flux_)

Full_TestData = storing(TestPLATEIDs,supers)
XTest,YTest,Test_z, Test_mag,TestAnd, TestIn, Testwavst, TestID = MLAData(Full_TestData,TestBinInfos,TestFlux, Testlog_wavst,TestANDMASK,TestINVAR)
plate_no = 0
all_Testwav=[]

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
XTest_plate=[]
print('Re-binning')
while plate_no<len(XTest):
    plateX = XTest[plate_no]
    plateMask = TestAnd[plate_no]
    plateInv = TestIn[plate_no]
    wave = all_Testwav[plate_no]
    r,w, r_weight = StandardRebin(plateX,wave,plateMask,plateInv ,Bin_Size )
    XTest_plate.append(r)
    plate_no=plate_no+1
    

p = 0

while p < len(XTest_plate):
    CurrentplateX = XTest_plate[p]
    CurrentplateY=YTest[p]
    n=0
    while n<len(CurrentplateX):
        X_Test.append(CurrentplateX[n][:800])
        Y_Test.append(CurrentplateY[n])
        n=n+1
    p=p+1













    
hiddenlayer_format = (13)
backprop_method = 'adam'
lr=0.0001
act =  'logistic'#'tanh'
mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
classes=[1,3,4,30]

object_total = 0
part_start = 0
increment = 1000
part_end = increment

while part_end<len(X_Full):
    try:
        X_current = X_Full[part_start:part_end]
        Y_current = Y_Full[part_start:part_end]
        mlp.partial_fit(X_current,Y_current,classes)
    except ValueError as V:
        del X_Full[part_start:part_end]
        del Y_Full[part_start:part_end]
        print(V)
        part_start=part_end
        print("Fitting: "+np.str((part_end/len(X_Full))*100))
        part_end+=increment
        object_total=part_end
        
    else:
        part_start=part_end
        print("Fitting: "+np.str((part_end/len(X_Full))*100))
        part_end+=increment
        object_total=part_end #+15 just in case of rounding
  

predictions = mlp.predict(np.array(X_Test))
star,star_starloc,star_lowzloc,star_galloc,star_highzloc = classification(1,Y_Test,predictions) 
lowz,lowz_starloc,lowz_loc,lowz_galloc,lowz_highzloc = classification(3,Y_Test,predictions)
gal,gal_starloc,gal_lowzloc,gal_galloc,gal_highzloc = classification(4,Y_Test,predictions)
highz,highz_starloc,highz_lowzloc,highz_galloc,highz_highzloc = classification(30,Y_Test,predictions)
File_Name = np.str(Bin_Size)#input("Please Enter File name: ")
d = open(File_Name+".txt", 'w')
#t1=["Files used",np.str(Spectra_Files), "\n"]
sp= "\n"
t2 = ["Bin Size = ",np.str(Bin_Size), "\n"]
t3 = ["Number of training objects = ",np.str(object_total), "\n"]
t4 = ["Number of testing objects = ",np.str(len(X_Test)), "\n"]
n1 = ["Structure of neural network: ", np.str(hiddenlayer_format),"\n"]
n2  = ["Backpropagation method used: ",np.str(backprop_method), "\n"]
n3  = ["Learning rate: ",np.str(lr), "\n"]
n4  = ["Activation Function: ",np.str(act), "\n"]
r1 = ["Results of Neural Network: ", "\n","\n"]
r2=["        ","       Star    Quasar  Galaxy  BAL ","\n",]
r3="Star           ",np.str(np.round(star[0]*100,2)),"%  ", np.str(np.round(star[1]*100,2)),"%  ",np.str(np.round(star[2]*100,2)),"%  ",np.str(np.round(star[3]*100,2)),"%","\n"
r4="Quasar z<2.1   ",np.str(np.round( lowz[0]*100,2)),"%  ", np.str(np.round( lowz[1]*100,2)),"%  ",np.str(np.round( lowz[2]*100,2)),"%  ",np.str(np.round( lowz[3]*100,2)),"%","\n"
r5="Galaxy         ",np.str(np.round( gal[0]*100,2)),"%  ", np.str(np.round( gal[1]*100,2)),"%  ",np.str(np.round( gal[2]*100,2)),"%  ",np.str(np.round( gal[3]*100,2)),"%","\n"
r6="Quasar z<2.1   ",np.str(np.round( highz[0]*100,2)),"%  ", np.str(np.round( highz[1]*100,2)),"%  ",np.str(np.round( highz[2]*100,2)),"%  ",np.str(np.round (highz[3]*100,2)),"%","\n"


#d.writelines(t1)
#d.writelines(sp)
#d.writelines(sp)
d.writelines(t2)
d.writelines(sp)
d.writelines(sp)
d.writelines(t3)
d.writelines(t4)
d.writelines(sp)
d.writelines(sp)
d.writelines(n1)
d.writelines(n2)
d.writelines(n3)
d.writelines(n4)
d.writelines(sp)
d.writelines(sp)
d.writelines(r1)
d.writelines(r2)
d.writelines(r3)
d.writelines(r4)
d.writelines(r5)
d.writelines(r6)

list.close()

d.close()


    
