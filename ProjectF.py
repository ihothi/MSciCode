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















def MLAData(Full_Data,BinInfos,Flux):
    
    Y = []
    X = []
    All_redshifts=[]
    All_Mag=[]
    plate_no = 0
    y=0
    while plate_no < len(Full_Data):
        #loading matched objects with sup, plate data
        CurrentSup_data = Full_Data[plate_no]
        CurrentBin = BinInfos[plate_no]
        CurrentFlux = Flux[plate_no]
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
                                else:
                                    Y.append(30)
                                    x_flux =(CurrentFlux[BinObj_No])
                                    x_flux=x_flux[:4600]
                                    X.append(x_flux)
                                    All_redshifts.append(CurrentSup.z)
                                    All_Mag.append(CurrentSup.Mag)
                                    
                            else: 
                                Y.append(CurrentSup.Class_p)
                                x_flux =(CurrentFlux[BinObj_No])
                                x_flux=x_flux[:4600]
                                X.append(x_flux)
                                All_redshifts.append(CurrentSup.z)
                                All_Mag.append(CurrentSup.Mag)
                BinObj_No=BinObj_No+1  
            Sup_obj=Sup_obj+1
        plate_no = plate_no+1

    return X,Y,All_redshifts,All_Mag


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