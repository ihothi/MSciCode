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
print("Retrieving Flux and Classification")
X,Y,Train_z, Train_mag,And, In, wavst, ID = MLAData(Full_Data,BinInfos,Flux, log_wavst,ANDMASK,INVAR)


X_Full = []
Y_Full = []
p = 0

while p < len(X):
    CurrentplateX = X[p]
    CurrentplateY=Y[p]
    n=0
    while n<len(CurrentplateX):
        if len(CurrentplateX[n])==0:
            n=n+1
        else:
            X_Full.append(CurrentplateX[n])
            Y_Full.append(CurrentplateY[n])
            n=n+1
    p=p+1
Train_No = len(X_Full)/4
X_Test=[]
Y_Test = []
i=0
while i< Train_No:
    a = randint(0, len(X_Full)-1)
    X_Test.append(X_Full[a])
    del X_Full[a]
    Y_Test.append(Y_Full[a])
    del Y_Full[a]
    i=i+1
    
    
hiddenlayer_format = (13)
backprop_method = 'sgd'
lr=0.0001
act =  'tanh' #'logistic'
mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About

#scaler = StandardScaler()  
#scaler.fit(X_Full)  
#X_Full = scaler.transform(X_Full)  
#X_Test = scaler.transform(X_Test)  


#mlp.fit(np.array(X_Full),np.array(Y_Full))
classes=[1,3,4,30]

object_total = 0
part_start = 0
increment = np.int((len(X_Full)/100))
part_end = increment

while object_total<len(X_Full):
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
        object_total=part_end+15
        
    else:
        part_start=part_end
        print("Fitting: "+np.str((part_end/len(X_Full))*100))
        part_end+=increment
        object_total=part_end+15 #+15 just in case of rounding
    
predictions = mlp.predict(np.array(X_Test))
star,star_starloc,star_lowzloc,star_galloc,star_highzloc = classification(1,Y_Test,predictions) 
lowz,lowz_starloc,lowz_loc,lowz_galloc,lowz_highzloc = classification(3,Y_Test,predictions)
gal,gal_starloc,gal_lowzloc,gal_galloc,gal_highzloc = classification(4,Y_Test,predictions)
highz,highz_starloc,highz_lowzloc,highz_galloc,highz_highzloc = classification(30,Y_Test,predictions)
File_Name = "Test_500plates"#input("Please Enter File name: ")
d = open(File_Name+".txt", 'w')
#t1=["Files used",np.str(Spectra_Files), "\n"]
sp= "\n"
#t2 = ["Files used to test: ",np.str(Spectra_TestFiles), "\n"]
t3 = ["Number of training objects = ",np.str(len(X_Full)), "\n"]
t4 = ["Number of testing objects = ",np.str(len(X_Test)), "\n"]
n1 = ["Structure of neural network: ", np.str(hiddenlayer_format),"\n"]
n2  = ["Backpropagation method used: ",np.str(backprop_method), "\n"]
n3  = ["Learning rate: ",np.str(lr), "\n"]
n4  = ["Activation Function: ",np.str(act), "\n"]
r1 = ["Results of Neural Network: ", "\n","\n"]
r2=["        ","       Star   z<2.1  Galaxy z<2.1 ","\n",]
r3="Star           ",np.str(np.round(star[0]*100,2)),"%  ", np.str(np.round(star[1]*100,2)),"%  ",np.str(np.round(star[2]*100,2)),"%  ",np.str(np.round(star[3]*100,2)),"%","\n"
r4="Quasar z<2.1   ",np.str(np.round( lowz[0]*100,2)),"%  ", np.str(np.round( lowz[1]*100,2)),"%  ",np.str(np.round( lowz[2]*100,2)),"%  ",np.str(np.round( lowz[3]*100,2)),"%","\n"
r5="Galaxy         ",np.str(np.round( gal[0]*100,2)),"%  ", np.str(np.round( gal[1]*100,2)),"%  ",np.str(np.round( gal[2]*100,2)),"%  ",np.str(np.round( gal[3]*100,2)),"%","\n"
r6="Quasar z<2.1   ",np.str(np.round( highz[0]*100,2)),"%  ", np.str(np.round( highz[1]*100,2)),"%  ",np.str(np.round( highz[2]*100,2)),"%  ",np.str(np.round (highz[3]*100,2)),"%","\n"


#d.writelines(t1)
#d.writelines(sp)
#d.writelines(t2)
#d.writelines(sp)
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

