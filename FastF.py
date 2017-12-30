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
        
        
def Faststoring(PLATEIDs,s,BinInfos,Flux,log_wavs,ANDMASK,INV, Bin_Size):
    All_Y=[]
    All_X = []
    All_redshifts=[]
    All_Mag=[]
    All_AND=[]
    All_Inv=[]
    All_Name=[]
    All_weight=[]
    Plate_Count =0
    supers = s
    bin_wav = Bin_Size*0.829026074968624
    while Plate_Count < len(PLATEIDs):
        Current_Plate = PLATEIDs[Plate_Count]
        Current_Bin = BinInfos[Plate_Count]
        Current_Flux =  Flux[Plate_Count]
        Current_log =  log_wavs[Plate_Count]
        Current_and =  ANDMASK[Plate_Count]
        Current_inv = INV[Plate_Count]
        CurrentObject = 0
        append_count=0
        wav_ratio = 10**0.0001
        cent_wav = 10**Current_log
        wavelength = []
        wavelength.append(cent_wav)
        current_wav = cent_wav
        while append_count < (4600-1):
            current_wav = current_wav*wav_ratio
            wavelength.append(current_wav)
            append_count=append_count+1
        platename_data = []
        currentSup_size =len(supers)
        while CurrentObject < currentSup_size:
            Current = supers[CurrentObject]
            no_match = True
            ##Seeing if current object is contained within a given plate
            if Current['PLATE'] == Current_Plate:
                 
                Object_ = Object(Current['SDSS_NAME'], Current['RA'], Current['Dec'], 
                         Current['Z_VI'], Current['CLASS_PERSON'],Current['PLATE'] ,
                         Current['MJD'], Current['FIBERID'],Current['PSFMAG'])
                ##We shall now see if a plates fits file contains this object
                BinObj_No = 0
                while BinObj_No<len(Current_Bin):
                    BinObj = Current_Bin[BinObj_No]
                    ObjAndM= Current_and[BinObj_No]
                    ObjInv = Current_inv[BinObj_No]
                    ##checking if the two match 
                    if no_match:
                        if BinObj['FIBERID'] == Object_.FiberID:
                            no_match = False
                            if Object_.Class_p == 3 or Object_.Class_p==30:
                                if Object_.z < 2.1:
                                    All_Y.append(1)
                                    flux =(Current_Flux[BinObj_No])
                                    flux=flux[:4600]
                                    rebin_flux=[]
                                    weight_=[]
                                    cwav=[]
                                    pix=0
                                    First = True
                                    bin_max = 3600
                                    while pix < len(x_flux):
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
                                                
                                            if pix< len(flux):
                                                w_ = 0
                                                m = Current_and[pix]
                                                wav = wavelength[pix]
                                                if m>0:
                                                    w_=0
                                                else:
                                                    w_ = Current_inv[pix]
                                                x_flux = x_flux +(flux[pix]*w_)
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
                                    All_weight.append(weight_)
                                    All_X.append(rebin_flux)
                                    All_redshifts.append(Object_.z)
                                    All_AND.append(ObjAndM)
                                    All_Inv.append(ObjInv)
                                    All_Mag.append(Object_.Mag)
                                    All_Name.append(Object_.name)
                                    
                                else:
                                    All_Y.append(3)
                                    flux =(Current_Flux[BinObj_No])
                                    flux=flux[:4600]
                                    rebin_flux=[]
                                    weight_=[]
                                    cwav=[]
                                    pix=0
                                    First = True
                                    bin_max = 3600
                                    while pix < len(x_flux):
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
                                                
                                            if pix< len(flux):
                                                w_ = 0
                                                m = Current_and[pix]
                                                wav = wavelength[pix]
                                                if m>0:
                                                    w_=0
                                                else:
                                                    w_ = Current_inv[pix]
                                                x_flux = x_flux +(flux[pix]*w_)
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
                                    All_weight.append(weight_)
                                    All_X.append(rebin_flux)
                                    All_redshifts.append(Object_.z)
                                    All_AND.append(ObjAndM)
                                    All_Inv.append(ObjInv)
                                    All_Mag.append(Object_.Mag)
                                    All_Name.append(Object_.name)
                        
                                    
                            else: 
                                if Object_.Class_p ==1:
                                    
                                    All_Y.append(0)
                                    flux =(Current_Flux[BinObj_No])
                                    flux=flux[:4600]
                                    rebin_flux=[]
                                    weight_=[]
                                    cwav=[]
                                    pix=0
                                    First = True
                                    bin_max = 3600
                                    while pix < len(x_flux):
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
                                                
                                            if pix< len(flux):
                                                w_ = 0
                                                m = Current_and[pix]
                                                wav = wavelength[pix]
                                                if m>0:
                                                    w_=0
                                                else:
                                                    w_ = Current_inv[pix]
                                                x_flux = x_flux +(flux[pix]*w_)
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
                                    All_weight.append(weight_)
                                    All_X.append(rebin_flux)
                                    All_redshifts.append(Object_.z)
                                    All_AND.append(ObjAndM)
                                    All_Inv.append(ObjInv)
                                    All_Mag.append(Object_.Mag)
                                    All_Name.append(Object_.name)
                                else:
                                    All_Y.append(2)
                                    flux =(Current_Flux[BinObj_No])
                                    flux=flux[:4600]
                                    rebin_flux=[]
                                    weight_=[]
                                    cwav=[]
                                    pix=0
                                    First = True
                                    bin_max = 3600
                                    while pix < len(x_flux):
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
                                                
                                            if pix< len(flux):
                                                w_ = 0
                                                m = Current_and[pix]
                                                wav = wavelength[pix]
                                                if m>0:
                                                    w_=0
                                                else:
                                                    w_ = Current_inv[pix]
                                                x_flux = x_flux +(flux[pix]*w_)
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
                                    All_weight.append(weight_)
                                    All_X.append(rebin_flux)
                                    All_redshifts.append(Object_.z)
                                    All_AND.append(ObjAndM)
                                    All_Inv.append(ObjInv)
                                    All_Mag.append(Object_.Mag)
                                    All_Name.append(Object_.name)
                BinObj_No=BinObj_No+1
                ## Once matched, no point
                ##As we have deleted an element 

            CurrentObject=CurrentObject+1 
        Plate_Count = Plate_Count + 1
        print("Plate Storage Completed: "+  np.str((Plate_Count*100)/len(PLATEIDs))+"%" )
    return All_Y, All_X, All_redshifts, All_Mag, All_AND,All_Inv,All_Name, All_weight 



