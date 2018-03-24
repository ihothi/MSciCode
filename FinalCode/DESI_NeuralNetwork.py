import numpy as np
from astropy.io import fits
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import os
from Project_Functions import MLAData,classification, Object,storing,MLADataBin, DESI
from random import randint

#Spectral Flux Bin size
Bin_Size = 10
#Rejection Threshhold
rejections = 500
slash =  os.path.normpath("/")

hiddenlayer_format = 80,80,80,80
backprop_method = 'lbfgs'#adam'
lr=0.0001
act = 'tanh'#  'logistic'



#Location of SDSS data
Platedir = os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"New_boss_data_rebinned_x10x"+slash)
plate_name = os.listdir(Platedir)




#Location of DESI data
desidir =  os.path.normpath(slash+"share"+slash+"splinter"+slash+"hothi"+slash+"Newdesi_sim_data")

#Opening and storing DESI data
desi_name = os.listdir(desidir)
desifile="FullDesi_Name.txt"
p = open(desifile, 'w')
for i in desi_name:
    p.write(i +'\n')
p.close()
DesiFiles=[]
with open(desifile) as f:
    DesiFiles = f.read().splitlines()
Desi_folder=[]
i=0
#len(DesiFiles)
while i<len(DesiFiles) :
    Desi_folder.append(DesiFiles[i])
    i=i+1
DesiFiber = []
FluxB = []
FluxR = []
FluxZ = []
WAVB = []
WAVR = []
WAVZ = []
    
MJDs = []
ANDMASKB=[]
INVARB=[]
ANDMASKZ=[]
INVARZ=[]
ANDMASKR=[]
INVARR=[]
Desi_table=[]
no_obj = []
for f in Desi_folder:
    file_list = os.listdir(desidir+slash+f)
    for l in file_list:
        c=desidir+slash+f+slash+l
        desiplate_ = fits.open(c,memmap=True)
        primhdu_desi = desiplate_[0]
        hdudesi =  desiplate_[1].data
        Desi_table.append(hdudesi)
        DesiFiber.append(primhdu_desi.header['FIBER'])
        no_obj.append(primhdu_desi.header['OBJNO'])
Bin_Size=10
i = 0
DesiIDs = []
Desi_Class = []
Desi_Z = []

while i < len(DesiFiber):
    plate_hdu  = Desi_table[i]
    obj = 0
    p_WAVB =[]
    p_WAVR=[]
    p_WAVZ=[]
    p_DesiIDs=[]
    p_FluxB=[]
    p_FluxR=[]
    p_FluxZ=[]
    p_ANDMASKB=[]
    p_INVARB=[]
    p_ANDMASKZ=[]
    p_INVARZ=[]
    p_ANDMASKR=[]
    p_INVARR=[]
    p_Desi_Class=[]
    p_DesiIDs=[]
    p_Desi_Z=[]
    print(len(plate_hdu))
    while obj <len(plate_hdu):
        currentobj = plate_hdu[obj]
        CurrentFluxB = currentobj[0]
        CurrentANDMASKB= currentobj[1]
        CurrentINVARB= currentobj[2]
        CurrentFluxR = currentobj[4]
        CurrentANDMASKR= currentobj[5]
        CurrentINVARR= currentobj[6]
        CurrentFluxZ = currentobj[8]
        CurrentANDMASKZ= currentobj[9]
        CurrentINVARZ= currentobj[10]
        CurrentClass = currentobj[12]
        CurrentZ = currentobj[13]
        CurrentID = currentobj[14]
        p_DesiIDs.append(CurrentID)
        p_FluxB.append(CurrentFluxB)
        p_FluxR.append(CurrentFluxR)
        p_FluxZ.append(CurrentFluxZ)
        p_ANDMASKB.append(CurrentANDMASKB)
        p_INVARB.append(CurrentINVARB)
        p_ANDMASKZ.append(CurrentANDMASKZ)
        p_INVARZ.append(CurrentINVARZ)
        p_ANDMASKR.append(CurrentANDMASKR)
        p_INVARR.append(CurrentINVARR)
        p_Desi_Class.append(CurrentClass)
        p_Desi_Z.append(CurrentZ)
        obj = obj +1
    WAVB.append(plate_hdu.field(3))
    WAVR.append(plate_hdu.field(7))
    WAVZ.append(plate_hdu.field(11))
    DesiIDs.append(p_DesiIDs)
    FluxB.append(p_FluxB)
    FluxR.append(p_FluxR)
    FluxZ.append(p_FluxZ)
    ANDMASKB.append(p_ANDMASKB)
    INVARB.append(p_INVARB)
    ANDMASKZ.append(p_ANDMASKZ)
    INVARZ.append(p_INVARZ)
    ANDMASKR.append(p_ANDMASKR)
    INVARR.append(p_INVARR)
    Desi_Class.append(p_Desi_Class)
    p_Desi_Z.append(p_Desi_Z)
    i = i+1




file="FullPlate_Name.txt"
p = open(file, 'w')
for i in plate_name:
    p.write(i +'\n')

p.close()


with open(file) as f:
    Files = f.read().splitlines() 
    
No_TrainPlates = 2200 ##Note there are only 1200 plates total      
i=0
Spectra_Files=[]
while i<2200:
    Spectra_Files.append(Files[i])
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


i=0    
x_f=[]
y_f=[]
DesiX_Full = []
DesiY_Full = []
minimumdesi = 100000 

#Binning DESI data   
while i < len(DesiIDs):
    print("Plate: ", i)
    plate_wavB=WAVB[i][:2380]
    plate_wavR=WAVR[i][:2116]
    plate_wavZ=WAVZ[i]
    plate_fluxB = FluxB[i]
    plate_fluxR = FluxR[i]
    plate_fluxZ = FluxZ[i]
    plate_andB = ANDMASKB[i]
    plate_andR = ANDMASKR[i]
    plate_andZ = ANDMASKZ[i]
    plate_inB = INVARB[i]
    plate_inR = INVARR[i]
    plate_inZ = INVARZ[i]
    plate_class = Desi_Class[i]
    k = no_obj[i]
    xfp,yfp,w,minim = DESI(plate_wavB,plate_wavR,plate_wavZ,plate_fluxB,plate_fluxR,plate_fluxZ,plate_andB,plate_andR,plate_andZ,plate_inB,plate_inR,plate_inZ,plate_class,Bin_Size,k,rejections)
    DesiX_Full.append(xfp)
    DesiY_Full.append(yfp)
    if minim < minimum:
        minimum = minim
    i+=1    
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
        plate_x.append(Currentplate_x[:minimum])
        plate_y.append(Currentplate_y)
        obj = obj +1
    i = i+1
        


i=2200
No_TestPlates = 470 ##Do not need many, say, <50
Spectra_TestFiles=[]
while i<2400 :
    Spectra_TestFiles.append(Files[i])
    i=i+1

#Opening and Storing SDSS data
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
        plate_xt.append(Currentplate_x[:minimum])
        plate_yt.append(Currentplate_y)
        obj = obj +1
    i = i+1
    
    




#Appending DESI flux to SDSS'   
for df in DesiX_Full:
    for x in df:
        plate_x.append(x[:minimum])
for dc in DesiY_Full:
    for y in dc:
        plate_y.append(y)


#Training algorithm 

mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
all_class=[0,1,2,3]

scaler = StandardScaler()  
scaler.fit(plate_x)  
plate_x = scaler.transform(plate_x)  
mlp.fit(plate_x, plate_y)
plate_xt = scaler.transform(plate_xt)  

#Predicting
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


#print(t1)
#print(sp)
#print(sp)
print(t2)
print(sp)
print(sp)
print(t3)
print(t4)
print(sp)
print(sp)
print(n1)
print(n2)
print(n3)
print(n4)
print(sp)
print(sp)
print(r1)
print(r2)
print(r3)
print(r4)
print(r5)
print(r6)
d.close()


    
