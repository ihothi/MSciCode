import numpy as np


class Object:
    
    def __init__():[]
   
    def __init__(self, SDSS_NAME, R, D, Z_VI, CLASS_P, p, mjd, fid,PSFMAG):
        #leaving out redshift for now
        R = round(R,2)
        D =  round(D,2)
        self.name = SDSS_NAME
        self.RA = R
        self.Dec = D
        self.z = Z_VI
        self.Class_p = CLASS_P
        self.Plate = p
        self.MJD = mjd
        self.FiberID = fid
        self.Mag = PSFMAG



class FluxStore:
    
    def __init__():[]
   
    def __init__(self,PF):
        #leaving out redshift for now
        self.PlateFlux = PF
        
        
        
        
def Rebin(Plate_hdu, Bin_Size ):
    
    rebin = []
    rebin_weight=[]
    obj=0
    while obj < len(Plate_hdu):
        current_flux = Plate_hdu[obj][0]
        current_mask = Plate_hdu[obj][2]
        current_inv = Plate_hdu[obj][3]
        rebin_flux=[]
        weight_=[]
        pix=0
        while pix < len(current_flux):
            print('yes')
            bin_no=0
            x_flux=0
            W=0
            while bin_no<Bin_Size:
                pix = pix+bin_no
                if pix< len(current_flux):
                    w_ = 0
                    m = current_mask[pix]
                    if m==0:
                        w_=0
                    else:
                        w_ = current_inv[pix]
                    x_flux = x_flux +(current_flux[pix]*w_)
                    W=W+w_
                    bin_no=bin_no+1
                else: 
                    bin_no=bin_no+1
            if W == 0:
                x_flux=0
            else: 
                x_flux = x_flux/W
            rebin_flux.append(x_flux)

            weight_.append(W)
        rebin.append(rebin_flux)
        rebin_weight.append(weight_)
        obj=obj+1
    return rebin, rebin_weight











def MLAData(Full_Data,BinInfos,Flux,log_wavs,ANDMASK, INV):
    
    
    All_Y=[]
    All_X = []
    All_redshifts=[]
    All_Mag=[]
    All_AND=[]
    All_Inv=[]
    All_Name=[]
    wav_logs=[]
    plate_no = 0
    y=0
    while plate_no < len(Full_Data):
        #loading matched objects with sup, plate data
        CurrentSup_data = Full_Data[plate_no]
        CurrentBin = BinInfos[plate_no]
        CurrentFlux = Flux[plate_no]
        CurrentAndM= ANDMASK[plate_no]
        CurrentInv = INV[plate_no]
        wav= log_wavs[plate_no]
        Plate_Y = []
        Plate_X = []
        Plate_AND=[]
        Plate_Inv=[]
        Plate_redshifts=[]
        Plate_Mag=[]
        Plate_Name=[]
        
        #first object is zeroth element
        Sup_obj =0 
        while Sup_obj < len(CurrentSup_data):
            #to stop double counting
            no_match = True
            #looking at an object that has been matched in sup list
            CurrentSup = CurrentSup_data[Sup_obj]            
            #going through each object in the bin
            BinObj_No = 0
            while BinObj_No<len(CurrentBin):
                BinObj = CurrentBin[BinObj_No]
                ObjAndM= CurrentAndM[plate_no]
                ObjInv = CurrentInv[plate_no]
                
                ##checking if the two match 
                if BinObj['FIBERID'] == CurrentSup.FiberID:
                    y=y+1 
                    if no_match:
                        if CurrentSup.Class_p == 0:
                            a=0
                        else:
    
                            no_match = False
                            if CurrentSup.Class_p == 3 or CurrentSup.Class_p==30:
                                if CurrentSup.z < 2.1:
                                    Plate_Y.append(3)
                                    x_flux =(CurrentFlux[BinObj_No])
                                    x_flux=x_flux[:4600]
                                    Plate_X.append(x_flux)
                                    Plate_redshifts.append(CurrentSup.z)
                                    Plate_AND.append(ObjAndM)
                                    Plate_Inv.append(ObjInv)
                                    Plate_Mag.append(CurrentSup.Mag)
                                    Plate_Name.append(CurrentSup.name)
                                    
                                else:
                                    Plate_Y.append(30)
                                    x_flux =(CurrentFlux[BinObj_No])
                                    x_flux=x_flux[:4600]
                                    Plate_X.append(x_flux)
                                    Plate_redshifts.append(CurrentSup.z)
                                    Plate_AND.append(ObjAndM)
                                    Plate_Inv.append(ObjInv)
                                    Plate_Mag.append(CurrentSup.Mag)
                                    Plate_Name.append(CurrentSup.name)
                        
                                    
                            else: 
                                Plate_Y.append(CurrentSup.Class_p)
                                x_flux =(CurrentFlux[BinObj_No])
                                x_flux=x_flux[:4600]
                                Plate_X.append(x_flux)
                                Plate_redshifts.append(CurrentSup.z)
                                Plate_AND.append(ObjAndM)
                                Plate_Inv.append(ObjInv)
                                Plate_Mag.append(CurrentSup.Mag)
                                Plate_Name.append(CurrentSup.name)
            
                BinObj_No=BinObj_No+1  
            Sup_obj=Sup_obj+1
        All_Y.append(Plate_Y)
        All_X.append(Plate_X)
        All_Name.append(Plate_Name)
        plate_no = plate_no+1
        wav_logs.append(wav)
        All_redshifts.append(Plate_redshifts)
        All_Mag.append(Plate_Mag)
        All_AND.append(Plate_AND)
        All_Inv.append(Plate_Inv)
    

    return All_X,All_Y,All_redshifts,All_Mag,All_AND,All_Inv,wav_logs,All_Name


def classification(objectclass, Trainingclass, prediction):
    ## Star, Quasar, Galaxy, BAL
    classi = []
    ##Location of object that is predicted to be a given classification
    loc = [0,0,0,0]
    star=0;
    #quasar w/ redshift <2.1
    qso=0;
    gal=0;
    #quasar w/ redshift >2.1
    bal=0;
    starloc=[];
    qsoloc=[];
    galloc=[];
    balloc=[];
    i=0
    while i< len(Trainingclass):
        currentobject = Trainingclass[i]
        currentpred = prediction[i]
        if objectclass==currentobject:
            if currentpred==1: 
                star=star+1
                starloc.append(i)
            elif currentpred==3: 
                qso=qso+1
                qsoloc.append(i)
            elif currentpred==4: 
                gal=gal+1
                galloc.append(i)
            elif currentpred==30: 
                bal=bal+1
                balloc.append(i)
        
        i=i+1
    
    classi.append(star/(star+qso+gal+bal))
    classi.append(qso/(star+qso+gal+bal))
    classi.append(gal/(star+qso+gal+bal))
    classi.append(bal/(star+qso+gal+bal))
    return classi,starloc,qsoloc,galloc,balloc


def storing(PLATEIDs,supers):
    CurrentPlate = 0
    Full_Data=[]
    Plate_Count =0
    while Plate_Count < len(PLATEIDs):
        Current_Plate = PLATEIDs[Plate_Count]
        CurrentObject = 0
        platename_data = []
        while CurrentObject < len(supers):
            Current = supers[CurrentObject]
            if Current['PLATE'] == Current_Plate:
                Object_ = Object(Current['SDSS_NAME'], Current['RA'], Current['Dec'], 
                         Current['Z_VI'], Current['CLASS_PERSON'],Current['PLATE'] ,
                         Current['MJD'], Current['FIBERID'],Current['PSFMAG'])
                platename_data.append(Object_)
            CurrentObject=CurrentObject+1    
        Plate_Count = Plate_Count + 1
        Full_Data.append(platename_data)
    return Full_Data

def MLADataBin(Full_Data,BinInfos,Flux,log_wavs, ANDMASK, INV, BIN_Size):
    
    
    All_Y=[]
    All_X = []
    All_redshifts=[]
    All_Mag=[]
    wav_logs=[]
    All_AND=[]
    All_Inv=[]
    plate_no = 0
    y=0
    while plate_no < len(Full_Data):
        #loading matched objects with sup, plate data
        CurrentSup_data = Full_Data[plate_no]
        CurrentBin = BinInfos[plate_no]
        CurrentFlux = Flux[plate_no]
        CurrentAndM= ANDMASK[plate_no]
        CurrentInv = INV[plate_no]
        wav= log_wavs[plate_no]
        Plate_Y = []
        Plate_X = []
        Plate_AND=[]
        Plate_Inv=[]
        
        #first object is zeroth element
        Sup_obj =0 
        while Sup_obj < len(CurrentSup_data):
            #to stop double counting
            no_match = True
            #looking at an object that has been matched in sup list
            CurrentSup = CurrentSup_data[Sup_obj]            
            #going through each object in the bin
            BinObj_No = 0
            while BinObj_No<len(CurrentBin):
                X=[]
                BinObj = CurrentBin[BinObj_No]
                ObjAndM= CurrentAndM[plate_no]
                ObjInv = CurrentInv[plate_no]
                
                ##checking if the two match 
                if BinObj['FIBERID'] == CurrentSup.FiberID:
                    y=y+1 
                    if no_match:
                        if CurrentSup.Class_p == 0:
                            a=0
                        else:
    
                            no_match = False
                            if CurrentSup.Class_p == 3 or CurrentSup.Class_p==30:
                                if CurrentSup.z < 2.1:
                                    Plate_Y.append(3)
                                    x_flux =(CurrentFlux[BinObj_No])
                                    x_flux=x_flux[:4600]
                                    All_redshifts.append(CurrentSup.z)
                                    All_Mag.append(CurrentSup.Mag)
                                    wav_logs.append(wav)
                                    pixno=0
                                    while pixno<len(x_flux):
                                        bin_count =0
                                        x_= 0
                                        w=0
                                        while bin_count<BIN_Size:
                                            x_ = x_+x_flux[pixno+bin_count]*(ObjInv[pixno+bin_count])
                                            w=w+ObjInv[pixno+bin_count]
                                        x_=x_/w
                                        X.append(x_)
                                    Plate_X.append(X)
                                        
                                        
                                            
                                            
                                            
                                else:
                                    Plate_Y.append(30)
                                    x_flux =(CurrentFlux[BinObj_No])
                                    x_flux=x_flux[:4600]
                                    All_redshifts.append(CurrentSup.z)
                                    All_Mag.append(CurrentSup.Mag)
                                    wav_logs.append(wav)
                                    while pixno<len(x_flux):
                                        bin_count =0
                                        x_= 0
                                        w=0
                                        while bin_count<BIN_Size:
                                            x_ = x_+x_flux[pixno+bin_count]*(ObjInv[pixno+bin_count])
                                            w=w+ObjInv[pixno+bin_count]
                                        x_=x_/w
                                        X.append(x_)
                                    Plate_X.append(X)
                                    
                            else: 
                                Plate_Y.append(CurrentSup.Class_p)
                                x_flux =(CurrentFlux[BinObj_No])
                                x_flux=x_flux[:4600]
                                All_redshifts.append(CurrentSup.z)
                                All_Mag.append(CurrentSup.Mag)
                                wav_logs.append(wav)
                                while pixno<len(x_flux):
                                    bin_count =0
                                    x_= 0
                                    w=0
                                    while bin_count<BIN_Size:
                                        x_ = x_+x_flux[pixno+bin_count]*(ObjInv[pixno+bin_count])
                                        w=w+ObjInv[pixno+bin_count]
                                    x_=x_/w
                                    X.append(x_)
                                Plate_X.append(X)
                BinObj_No=BinObj_No+1  
                Plate_X
            Sup_obj=Sup_obj+1
        All_Y.append(Plate_Y)
        All_X.append(Plate_X)
        plate_no = plate_no+1

    return All_X,All_Y,All_redshifts,All_Mag,wav_logs


def MLADataTest(Full_Data,BinInfos,Flux,log_wavs):
    
    
    Plate_Y = []
    Plate_X = []
    All_redshifts=[]
    All_Mag=[]
    wav_logs=[]
    
    plate_no = 0
    y=0
    while plate_no < len(Full_Data):
        Y = []
        X = []
        #loading matched objects with sup, plate data
        CurrentSup_data = Full_Data[plate_no]
        CurrentBin = BinInfos[plate_no]
        CurrentFlux = Flux[plate_no]
        wav= log_wavs[plate_no]
        #first object is zeroth element
        Sup_obj =0 
        while Sup_obj < len(CurrentSup_data):
            #to stop double counting
            no_match = True
            #looking at an object that has been matched in sup list
            CurrentSup = CurrentSup_data[Sup_obj]            
            #going through each object in the bin
            BinObj_No = 0
            while BinObj_No<len(CurrentBin):
                BinObj = CurrentBin[BinObj_No]
                
                ##checking if the two match 
                if BinObj['FIBERID'] == CurrentSup.FiberID:
                    y=y+1 
                    if no_match:
                        if CurrentSup.Class_p == 0:
                            a=0
                        else:
    
                            no_match = False
                            if CurrentSup.Class_p == 3 or CurrentSup.Class_p==30:
                                if CurrentSup.z < 2.1:
                                    Y.append(3)
                                    x_flux =(CurrentFlux[BinObj_No])
                                    x_flux=x_flux[:4600]
                                    X.append(x_flux)
                                    All_redshifts.append(CurrentSup.z)
                                    All_Mag.append(CurrentSup.Mag)
                                    wav_logs.append(wav)
                                else:
                                    Y.append(30)
                                    x_flux =(CurrentFlux[BinObj_No])
                                    x_flux=x_flux[:4600]
                                    X.append(x_flux)
                                    All_redshifts.append(CurrentSup.z)
                                    All_Mag.append(CurrentSup.Mag)
                                    wav_logs.append(wav)
                                    
                            else: 
                                Y.append(CurrentSup.Class_p)
                                x_flux =(CurrentFlux[BinObj_No])
                                x_flux=x_flux[:4600]
                                X.append(x_flux)
                                All_redshifts.append(CurrentSup.z)
                                All_Mag.append(CurrentSup.Mag)
                                wav_logs.append(wav)
                BinObj_No=BinObj_No+1  
            Sup_obj=Sup_obj+1
        plate_no = plate_no+1
        Plate_X.append(X)
        Plate_Y.append(Y)

    return Plate_X,Plate_Y,All_redshifts,All_Mag,wav_logs

