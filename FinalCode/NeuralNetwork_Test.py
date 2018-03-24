import numpy as np
from astropy.io import fits
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
from Project_Functions import classification
from sklearn.externals import joblib

#Path to Training data
slash =  os.path.normpath("/")
Platedir = os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"New_boss_data_rebinned_x10x"+slash)
plate_name = os.listdir(Platedir)

#Path To model 
model = '/home/ihothi/optNN.pkl'
#path to scaler
scale_load = '/home/ihothi/scaler.save' 


#Results Storing File 
File_Name = 'NN_Results'



file="FullPlate_Name.txt"
p = open(file, 'w')
for i in plate_name:
    p.write(i +'\n')
p.close()
Files = []
with open(file) as f:
    Files = f.read().splitlines() 
    

##Do not need many, say, <50
Spectra_TestFiles=[]
for f in Files :
    Spectra_TestFiles.append(f)

#Opening testing plate files
minimum=0
wav_testlog=[]
Full_Testtable = []
TestPLATEIDs = []
for f in Spectra_TestFiles:
    file_list = os.listdir(Platedir+slash+f)
    for l in file_list:
        c=Platedir+slash+f+slash+l
        plate_ = fits.open(c,memmap=True)
        primhdu_ = plate_[0]
        hdu =  plate_[1].data
        Full_Testtable.append(hdu)
        TestPLATEIDs.append(primhdu_.header['Plate'])
        wav_testlog.append(primhdu_.header['LogWav'])
        minimum = primhdu_.header['Min']
        
plate_rt  = []
plate_nt  = []
plate_ct  = []
plate_xt  = []
plate_yt  = []
i = 0

#Storing plate data
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
        plate_xt.append(Currentplate_x[:minimum])
        plate_yt.append(Currentplate_y)
        obj = obj +1
    i = i+1



#Loading model and scaler
mlp = joblib.load(model)
scaler = joblib.load(scale_load)

#Feature scaling
plate_xt = scaler.transform(plate_xt)  

#Predicting
predictions = mlp.predict(plate_xt)

star,star_starloc,star_lowzloc,star_galloc,star_highzloc = classification(0,plate_yt,predictions) 
lowz,lowz_starloc,lowz_loc,lowz_galloc,lowz_highzloc = classification(1,plate_yt,predictions)
gal,gal_starloc,gal_lowzloc,gal_galloc,gal_highzloc = classification(2,plate_yt,predictions)
highz,highz_starloc,highz_lowzloc,highz_galloc,highz_highzloc = classification(3,plate_yt,predictions)

d = open(File_Name+".txt", 'w')

r2=["        ","       Star   Quasar z<2.  Galaxy  Quasar z>2.1  ","\n",]
r3="Star           ",np.str(np.round(star[0]*100,2)),"%  ", np.str(np.round(star[1]*100,2)),"%  ",np.str(np.round(star[2]*100,2)),"%  ",np.str(np.round(star[3]*100,2)),"%","\n"
r4="Quasar z<2.1   ",np.str(np.round( lowz[0]*100,2)),"%  ", np.str(np.round( lowz[1]*100,2)),"%  ",np.str(np.round( lowz[2]*100,2)),"%  ",np.str(np.round( lowz[3]*100,2)),"%","\n"
r5="Galaxy         ",np.str(np.round( gal[0]*100,2)),"%  ", np.str(np.round( gal[1]*100,2)),"%  ",np.str(np.round( gal[2]*100,2)),"%  ",np.str(np.round( gal[3]*100,2)),"%","\n"
r6="Quasar z>2.1   ",np.str(np.round( highz[0]*100,2)),"%  ", np.str(np.round( highz[1]*100,2)),"%  ",np.str(np.round( highz[2]*100,2)),"%  ",np.str(np.round (highz[3]*100,2)),"%","\n"

d.writelines(r2)
d.writelines(r3)
d.writelines(r4)
d.writelines(r5)
d.writelines(r6)

list.close()

d.close()



    
