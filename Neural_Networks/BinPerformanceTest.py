import numpy as np
from astropy.io import fits
import os
from ProjectF import MLAData,classification, Object,storing
import matplotlib
import matplotlib.pyplot as plt
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from random import randint
import time

def StandardRebin(plateX, wavelength,ANDMASK,INVAR ,Bin_Size,minimum ):
    
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
            rebin_flux = rebin_flux[:800]
        
            

            weight_.append(W)
        rebinwav.append(cwav)
        rebin.append(rebin_flux)
        rebin_weight.append(weight_)
        if len(rebin_flux)<minimum:
            minimum=len(rebin_flux)
        obj=obj+1
    return rebin, rebin_weight,rebinwav, minimum


Bins = [1,3,10,30]

for Bin_Size in Bins:
    
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
        
    Spectra_Files=[]        
    i=500
    Spectra_Files=[]
    while i< 600:
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
    minimum=1000
    print('Re-binning Test')
    while plate_no<len(X):
        plateX = X[plate_no]
        plateMask = And[plate_no]
        plateInv = In[plate_no]
        wave = all_Testwav[plate_no]
        r,w, r_lambda,m = StandardRebin(plateX,wave,plateMask,plateInv ,Bin_Size,minimum )
        minimum=m
        Xbin_plate.append(r)
        w_plate.append(w)
        l_plate.append(r_lambda)
        plate_no=plate_no+1
    
            
    print("Fitting MLA") 
    
    hiddenlayer_format = 13
    backprop_method = 'lbfgs'#adam'
    lr=0.0001
    act =  'logistic'#'tanh'
    mlp = MLPClassifier(hidden_layer_sizes=hiddenlayer_format,max_iter=500, solver = backprop_method,learning_rate_init=lr,activation=act) ##Think About
    all_class=[0,1,2,3]
    
    
    
        #while part_end<len(X_Full):
        # try:
            #    X_current = X_Full[part_start:part_end]
            #     Y_current = Y_Full[part_start:part_end]
            #      mlp.partial_fit(X_current,Y_current,classes=all_class)
            #   except ValueError as V:
            #  del X_Full[part_start:part_end]
            #   del Y_Full[part_start:part_end]
            #    print(V)
   # else:
  #      part_start=part_end
 #       part_end+=increment
#        object_total=part_end #+15 just in case of rounding
  
    print("Predicting")
     
     
    Spectra_TestFiles=[]        
    i=800
    Spectra_TestFiles=[]
    while i< 850:
        Spectra_TestFiles.append(Files[i])
        i=i+1
    PLATEIDsTest = []
    BinInfosTest = []
    FluxTest = []
    MJDsTest = []
    log_wavstTest=[]
    ORMASKTest=[]
    ANDMASKTest=[]
    INVARTest=[]
    print("Opening Files")
    for f in Spectra_TestFiles:
        file_list = os.listdir(Platedir+slash+f)
        for l in file_list:
            if 'spPlate' in l and ".fits"in l: 
                c=Platedir+slash+f+slash+l
                plate_ = fits.open(c,memmap=True)
                Bin_info_ = plate_[5].data
                Flux_ = plate_[0].data
                primhdu_ = plate_[0]
                PLATEIDsTest.append(primhdu_.header['PLATEID'])
                ORMASKTest.append( plate_[3].data)
                ANDMASKTest.append( plate_[2].data)
                INVARTest.append( plate_[1].data)
                log_wavstTest.append(primhdu_.header['COEFF0'])
                MJDsTest.append(primhdu_.header['MJD'])
                BinInfosTest.append(Bin_info_)
                FluxTest.append(Flux_)
        

    print("Restoring Data")
    Full_DataTest = storing(PLATEIDsTest,supers)
    XTest,YTest,Train_zTest, Train_magTest,AndTest, InTest, wavstTest, IDTest, MJTest, MatchedPlatesTest = MLAData(Full_DataTest,BinInfosTest,FluxTest, log_wavstTest,ANDMASKTest,INVARTest,MJDsTest,PLATEIDsTest)
    plate_no = 0
    all_Testwav=[]
    wav_ratio = 10**0.0001
    
    while plate_no < len(wavstTest): 
        append_count=0
        cent_wav = 10**wavstTest[plate_no]
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
    Xbin_plateTest=[]
    w_plateTest=[]
    l_plateTest=[]
    print('Re-binning Test')
    while plate_no<len(XTest):
        plateXTest = XTest[plate_no]
        plateMaskTest = AndTest[plate_no]
        plateInvTest = InTest[plate_no]
        waveTest = all_Testwav[plate_no]
        rTest,wTest, r_lambdaTest,m = StandardRebin(plateXTest,waveTest,plateMaskTest,plateInvTest ,Bin_Size,minimum )
        minimum=m
        Xbin_plateTest.append(rTest)
        w_plateTest.append(wTest)
        l_plateTest.append(r_lambdaTest)
        plate_no=plate_no+1
        
    XFull = []
    YFull = []
    plate_count=0
    while plate_count < len(Xbin_plate):
        currentp_x = Xbin_plate[plate_count]
        platey = Y[plate_count]
        objcount=0
        while objcount <len(currentp_x):
            XFull.append(currentp_x[objcount][:minimum])
            YFull.append(platey[objcount])
            objcount=objcount+1
        plate_count=plate_count+1
    
    XFullTest = []
    YFullTest = []
    plate_count=0
    while plate_count <len(Xbin_plateTest):
        currentp_x = Xbin_plateTest[plate_count]
        platey = YTest[plate_count]
        objcount=0
        while objcount <len(currentp_x):
            XFullTest.append(currentp_x[objcount][:minimum])
            YFullTest.append(platey[objcount])
            objcount=objcount+1
        plate_count=plate_count+1
        
        
    scaler = StandardScaler()  
    scaler.fit(XFull)  
    XFull = scaler.transform(XFull)  
    
    print("Fitting")
    starttrain = time.time()
    mlp.fit(XFull, YFull)
    endtrain = time.time()
    print("time taken to train for :"+np.str(Bin_Size)+" bin size" + np.str(endtrain-starttrain)+" seconds" )
    XFullTest = scaler.transform(XFullTest) 
       
        
        
    teststart = time.time()
    predictions = mlp.predict(XFullTest)
    testend = time.time()
    print("time taken to test for :"+np.str(Bin_Size)+" bin size" +np.str(testend-teststart)+" seconds" )
    star,star_starloc,star_lowzloc,star_galloc,star_highzloc = classification(0,YFullTest,predictions) 
    lowz,lowz_starloc,lowz_loc,lowz_galloc,lowz_highzloc = classification(1,YFullTest,predictions)
    gal,gal_starloc,gal_lowzloc,gal_galloc,gal_highzloc = classification(2,YFullTest,predictions)
    highz,highz_starloc,highz_lowzloc,highz_galloc,highz_highzloc = classification(3,YFullTest,predictions)
    

        #t1=["Files used",np.str(Spectra_Files), "\n"]
    sp= "\n"
    t2 = ["Bin Size = "+np.str(Bin_Size)+ "\n"]
    t3 = ["Number of training objects = "+np.str(len(YFull))+ "\n"]
    t4 = ["Number of testing objects = "+np.str(len(YFullTest))+ "\n"]
    n1 = ["Structure of neural network: "+ np.str(hiddenlayer_format)+"\n"]
    n2  = ["Backpropagation method used: "+np.str(backprop_method)+ "\n"]
    n3  = ["Learning rate: "+np.str(lr)+ "\n"]
    n4  = ["Activation Function: "+np.str(act)+ "\n"]
    r1 = ["Results of Neural Network: "+ "\n"+"\n"]
    r2=["        "+"       Star    Quasar  Galaxy  BAL "+"\n"]
    r3="Star           "+np.str(np.round(star[0]*100,2))+"%  "+ np.str(np.round(star[1]*100,2)),"%  ",np.str(np.round(star[2]*100,2)),"%  ",np.str(np.round(star[3]*100,2)),"%","\n"
    r4="Quasar z<2.1   "+np.str(np.round( lowz[0]*100,2))+"%  "+ np.str(np.round( lowz[1]*100,2)),"%  ",np.str(np.round( lowz[2]*100,2)),"%  ",np.str(np.round( lowz[3]*100,2)),"%","\n"
    r5="Galaxy         "+np.str(np.round( gal[0]*100,2))+"%  "+ np.str(np.round( gal[1]*100,2)),"%  ",np.str(np.round( gal[2]*100,2)),"%  ",np.str(np.round( gal[3]*100,2)),"%","\n"
    r6="Quasar z<2.1   "+np.str(np.round( highz[0]*100,2))+"%  "+ np.str(np.round( highz[1]*100,2)),"%  ",np.str(np.round( highz[2]*100,2)),"%  ",np.str(np.round (highz[3]*100,2)),"%","\n"
    
    
    #prnt(t1)
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
            
        