import numpy as np

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
                            Y.append(CurrentSup.Class_p)
                            x_flux =(CurrentFlux[BinObj_No])
                            x_flux=x_flux[:4600]
                            X.append(x_flux)
                            All_redshifts.append(CurrentSup_test.z)
                            All_Mag.append(CurrentSup_test.Mag)
                BinObj_No=BinObj_No+1  
            Sup_obj=Sup_obj+1
        plate_no = plate_no+1

    return X,Y,All_redshifts,All_Mag


def classification(objectclass, Trainingclass, prediction):
    ## Star, Quasar, Galaxy, BAL
    classi = []
    star=0;
    qso=0;
    gal=0;
    bal=0;
    i=0
    while i< len(Trainingclass):
        currentobject = Trainingclass[i]
        if objectclass==1:
            if currentobject==1: 
                star=star+1
            elif currentobject==3: 
                qso=qso+1
            elif currentobject==4: 
                gal=gal+1
            elif currentobject==30: 
                bal=bal+1
                

        elif objectclass==3:
            if currentobject==1: 
                star=star+1
            elif currentobject==3: 
                qso=qso+1
            elif currentobject==4: 
                gal=gal+1
            elif currentobject==30: 
                bal=bal+1
                
        elif objectclass==4:
            if currentobject==1: 
                star=star+1
            elif currentobject==3: 
                qso=qso+1
            elif currentobject==4: 
                gal=gal+1
            elif currentobject==30: 
                bal=bal+1
                

        elif objectclass==30:
            if currentobject==1: 
                star=star+1
            elif currentobject==3: 
                qso=qso+1
            elif currentobject==4: 
                gal=gal+1
            elif currentobject==30: 
                bal=bal+1
        
        i=i+1
    
    classi.append(star)
    classi.append(qso)
    classi.append(gal)
    classi.append(bal)
    return classi
        
                
    
    

