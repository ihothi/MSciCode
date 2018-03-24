import numpy as np
from astropy.io import fits
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
from ProjectF import MLAData,classification, Object,storing,MLADataBin
from random import randint
from sklearn.externals import joblib
import time

# Hidden layer format (and number of nodes in each)
hiddenlayer_format = 80,80,80,80
#Back-propagation method
backprop_method = 'lbfgs'
#Learning Rate
lr=0.0001
#Activation Function
act = 'tanh'


#Enter path to data 
slash =  os.path.normpath("/")
Platedir = os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"New_boss_data_rebinned_x10x"+slash)
#List of plates
plate_name = os.listdir(Platedir)
file="NeuralNetwork_Plates.txt"
p = open(file, 'w')
for i in plate_name:
    p.write(i +'\n')
p.close()


with open(file) as f:
    Files = f.read().splitlines() 
    
No_TrainPlates = 2200 ##Note there are only 1200 plates total      
i=0
#Creating a list of training plates
Spectra_Files=[]
for f in Files:
    Spectra_Files.append(f)
    i=i+1

wav_log=[]
Full_table = []
PLATEIDs = []
minimum =0
#Opening the plate data
for f in Spectra_Files:
    file_list = os.listdir(Platedir+slash+f)
    for l in file_list:
        c=Platedir+slash+f+slash+l
        plate_ = fits.open(c,memmap=True)
        primhdu_ = plate_[0]
        hdu =  plate_[1].data
        Full_table.append(hdu)
        PLATEIDs.append(primhdu_.header['Plate'])
        wav_log.append(primhdu_.header['LogWav'])
        minimum = primhdu_.header['Min']
        
plate_r  = []
plate_n  = []
plate_c  = []
plate_x  = []
plate_y  = []
i = 0

#Storing the plate data
while i < len(PLATEIDs):
    plate_hdu  = Full_table[i]
    obj = 0
    while obj <len(plate_hdu):
        currentobj = plate_hdu[obj]
        Currentplate_z = currentobj[4]
        Currentplate_y = currentobj[1]
        Currentplate_x = currentobj[0]
        Currentplate_n = currentobj[5]
        plate_r.append(Currentplate_z)
        plate_n.append(Currentplate_n)
        plate_x.append(Currentplate_x[:minimum])
        plate_y.append(Currentplate_y)
        obj = obj +1
    i = i+1
        
hiddenlayer_format = 10
mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
all_class=[0,1,2,3]

scaler = StandardScaler()  
scaler.fit(plate_x)  
#Saving Feature Scaling
plate_x = scaler.transform(plate_x)  
#Fitting model

st = time.time()
mlp.fit(plate_x, plate_y)
et = time.time()
#Saving Model
print("Time: ", (et-st))


hiddenlayer_format = 20
mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
all_class=[0,1,2,3]

scaler = StandardScaler()  
scaler.fit(plate_x)  
#Saving Feature Scaling
plate_x = scaler.transform(plate_x)  
#Fitting model

st = time.time()
mlp.fit(plate_x, plate_y)
et = time.time()
#Saving Model
print("Time: ", (et-st))


hiddenlayer_format = 40
mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
all_class=[0,1,2,3]

scaler = StandardScaler()  
scaler.fit(plate_x)  
#Saving Feature Scaling
plate_x = scaler.transform(plate_x)  
#Fitting model

st = time.time()
mlp.fit(plate_x, plate_y)
et = time.time()
#Saving Model
print("Time: ", (et-st))


hiddenlayer_format = 80
mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
all_class=[0,1,2,3]

scaler = StandardScaler()  
scaler.fit(plate_x)  
#Saving Feature Scaling
plate_x = scaler.transform(plate_x)  
#Fitting model

st = time.time()
mlp.fit(plate_x, plate_y)
et = time.time()
#Saving Model
print("Time: ", (et-st))