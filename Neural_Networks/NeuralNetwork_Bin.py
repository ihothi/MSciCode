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
    Files = f.read().splitlines() 
    
No_TrainPlates = 1150 ##Note there are only 1200 plates total      
i=0
Spectra_Files=[]
while i<No_TrainPlates :
    a = randint(0, 1200)
    
    if a == 64:
        ##Dummy
        b=0
    else:
        Spectra_Files.append(Files[a])
        i=i+1

wav_log=[]
Full_table = []
PLATEIDs = []
minimum =0
print("Opening Files")
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
        plate_x.append(Currentplate_x[:Min])
        plate_y.append(Currentplate_y)
        obj = obj +1
    i = i+1
        


i=0
No_TestPlates = 20 ##Do not need many, say, <50
Spectra_TestFiles=[]
while i<No_TestPlates :
    a = randint(0, 1200)
    
    if a == 64:
        ##Dummy
        b=0
    else:
        Spectra_TestFiles.append(Files[a])
        i=i+1


wav_testlog=[]
Full_Testtable = []
TestPLATEIDs = []
print("Opening Files")
for f in Spectra_TestFiles:
    file_list = os.listdir(Platedir+slash+f)
    for l in file_list:
        c=Platedir+slash+f+slash+l
        print(c)
        plate_ = fits.open(c,memmap=True)
        primhdu_ = plate_[0]
        hdu =  plate_[1].data
        Full_Testtable.append(hdu)
        TestPLATEIDs.append(primhdu_.header['Plate'])
        wav_testlog.append(primhdu_.header['LogWav'])
        
plate_rt  = []
plate_nt  = []
plate_ct  = []
plate_xt  = []
plate_yt  = []
i = 0
while i < len(TestPLATEIDs):
    plate_hdu  = Full_Testtable[i]
    obj = 0
    while obj <len(plate_hdu):
        currentobj = plate_hdu[obj]
        Currentplate_z = currentobj[4]
        Currentplate_y = currentobj[1]
        Currentplate_x = currentobj[0]
        Currentplate_n = currentobj[5]
        plate_rt.append(Currentplate_z)
        plate_nt.append(Currentplate_n)
        plate_xt.append(Currentplate_x[:800])
        plate_yt.append(Currentplate_y)
        obj = obj +1
    i = i+1




print("Fitting MLA")    
hiddenlayer_format = 13
backprop_method = 'lbfgs'#adam'
lr=0.0001
act =  'logistic'#'tanh'
mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
all_class=[0,1,2,3]

scaler = StandardScaler()  
scaler.fit(plate_x)  
plate_x = scaler.transform(plate_x)  

print("Fitting")
mlp.fit(plate_x, plate_y)
     object_total=part_end #+15 just in case of rounding
  
print("Predicting")
plate_xt = scaler.transform(plate_xt)  

predictions = mlp.predict(plate_xt)
star,star_starloc,star_lowzloc,star_galloc,star_highzloc = classification(0,plate_yt,predictions) 
lowz,lowz_starloc,lowz_loc,lowz_galloc,lowz_highzloc = classification(1,plate_yt,predictions)
gal,gal_starloc,gal_lowzloc,gal_galloc,gal_highzloc = classification(2,plate_yt,predictions)
highz,highz_starloc,highz_lowzloc,highz_galloc,highz_highzloc = classification(3,plate_yt,predictions)
File_Name = np.str(Bin_Size)+"_1000plates"#input("Please Enter File name: ")
d = open(File_Name+".txt", 'w')
#t1=["Files used",np.str(Spectra_Files), "\n"]
sp= "\n"
t2 = ["Bin Size = ",np.str(Bin_Size), "\n"]
t3 = ["Number of training objects = ",np.str(len(plate_y)), "\n"]
t4 = ["Number of testing objects = ",np.str(len(plate_yt)), "\n"]
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



d.close()


    
