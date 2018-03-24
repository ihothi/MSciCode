import numpy as np

# Creating an object to store a Superset objects' data
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

def StandardRebin(plateX, wavelength,ANDMASK,INVAR ,Bin_Size,reject_no,plateY,z0,pid):  
    ''' Input: Objects Fluxes, Wavelength grid, ANDMASK, INVAR, size of bin, 
               Pixel rejection threshold, classifications, redshifts, plate ID  
    
        Output: Bin flux, INVAR of bin, wavelength of bin centre, 
                Minimum number of bins, Object class,  
        '''
   
    # Defining the return values
    rebin = []
    rebin_weight=[]
    rebin_class = []
    rebin_z = []
    rebin_id = []
    
    #Finding the minimum number of bins
    minimum=1000000
    # Ratio between bin size and starting-wavelength of bin
    bin_wav = Bin_Size*0.00023
    # For the binned wavelengths
    rebinwav=[]
    
    obj=0
    while obj < len(plateX):
        # initialising the number of pixels rejected variable 
        rej=0
        #Initialising object parameters
        current_flux = plateX[obj]
        Current_class = plateY[obj]
        current_z = z0[obj]
        current_id = pid[obj]
        current_mask = ANDMASK[obj]
        current_inv = INVAR[obj]
        rebin_flux=[]
        weight_=[]
        cwav=[]
        pix=0
        First = True
        bin_max = 3600
        
        # Going through each pixel in the spectrum
        while pix < len(current_flux):
            x_flux=0
            W=0
            wav=0
            # If the first wavelength, we need to first add its bin centre
            if First: 
                cwav.append(bin_max+((bin_max*bin_wav)*0.5))
                bin_max = bin_max+(bin_max*bin_wav)
                
                First = False
            # For the rest of the pixels
            else: 
                bin_max = bin_max+(bin_max*bin_wav)
                cwav.append(bin_max+((bin_max*bin_wav)*0.5))
            #Binning all the fluxes of wavelengths in a given bin
            while wav<bin_max:
                
                if pix< len(current_flux):
                    w_ = 0
                    m = current_mask[pix]
                    wav = wavelength[pix]
                    #rejecting pixels with ANDMASK = 0
                    if m>0:
                        rej+=1
                        w_=0
                    else:
                        w_ = current_inv[pix]
                    x_flux = x_flux +(current_flux[pix]*w_)
                    W=W+w_
                    pix=pix+1
                # To stop errors, if all the pixels in a spectrum have been
                # iterated over, exit loop
                else: 
                    pix=pix+1
                    wav=bin_max
                    
            # If the weight is zero, set flux of bin to zero
            if W == 0:
                x_flux=0
            
            else: 
                x_flux = x_flux/W
            rebin_flux.append(x_flux)
            weight_.append(W)
        #Checking If the number of rejected pixels is below the defined 
        #threshold
        if rej <reject_no:
            rebinwav.append(cwav)
            rebin.append(rebin_flux)
            rebin_weight.append(weight_)
            rebin_class.append(Current_class)
            rebin_z.append(current_z)
            rebin_id.append(current_id)
            #Finding the minimum number of bins
            if len(rebin_flux)<minimum:
                minimum=len(rebin_flux)
        obj=obj+1
    return rebin, rebin_weight,rebinwav,minimum,rebin_class,rebin_z,rebin_id





def MLAData(Full_Data,BinInfos,Flux,log_wavs,ANDMASK, INV, MJDs, PlateIDs):
    
    ''' Input: Superset data, Plate object data, logarithm of wavelengths, ANDMASK, MJDs, Plate IDs 
    
        Output: object fluxs, object classifications, object redshifts, 
                object magnitude, object ANDMASKS,object INVAR,
                object logarithm of wavelength, object name,
                object MJDs, PlateID
        '''
    # Initialising output variables
    All_Y=[]
    All_X = []
    All_redshifts=[]
    All_Mag=[]
    All_AND=[]
    All_Inv=[]
    All_Name=[]
    All_MJDs = []
    MatchedPlates = []
    wav_logs=[]
    plate_no = 0
    while plate_no < len(Full_Data):
        #loading matched objects with sup, plate data
        CurrentSup_data = Full_Data[plate_no]
        CurrentBin = BinInfos[plate_no]
        CurrentFlux = Flux[plate_no]
        CurrentAndM= ANDMASK[plate_no]
        CurrentInv = INV[plate_no]
        wav= log_wavs[plate_no]
        currentmjd = MJDs[plate_no]
        Plate_Y = []
        Plate_X = []
        Plate_AND=[]
        Plate_Inv=[]
        Plate_redshifts=[]
        Plate_Mag=[]
        Plate_Name=[]
        pid=PlateIDs[plate_no]

        plate_match = False
        
        #first object is zeroth element
        Sup_obj =0 
        #Going through superset objects
        while Sup_obj < len(CurrentSup_data):
            #to stop double counting
            no_match = True
            #looking at an object that has been matched in sup list
            CurrentSup = CurrentSup_data[Sup_obj]            
            #going through each object in the bin
            BinObj_No = 0
            while BinObj_No<len(CurrentBin):
                BinObj = CurrentBin[BinObj_No]
                ObjAndM= CurrentAndM[BinObj_No]
                ObjInv = CurrentInv[BinObj_No]
                
                ##checking if the two match 
                if BinObj['FIBERID'] == CurrentSup.FiberID:
                    if currentmjd==CurrentSup.MJD :
                        y+=1

                        plate_match = True 
                        if no_match:
                            #This class is assigned to objects that have not been inspected
                            # So we shall ignore
                            if CurrentSup.Class_p == 0:
                                a=0
                            else:
                                no_match = False
                                # Checking if it is a quasar to divide into
                                # redshift based classifications
                                if CurrentSup.Class_p == 3 or CurrentSup.Class_p==30:
                                    if CurrentSup.z < 2.1:
                                        Plate_Y.append(1)
                                        x_flux =(CurrentFlux[BinObj_No])
                                        x_flux=x_flux
                                        Plate_X.append(x_flux)
                                        Plate_redshifts.append(CurrentSup.z)
                                        Plate_AND.append(ObjAndM)
                                        Plate_Inv.append(ObjInv)
                                        Plate_Mag.append(CurrentSup.Mag)
                                        Plate_Name.append(CurrentSup.name)
                                       
                                        
                                    else:
                                        Plate_Y.append(3)
                                        x_flux =(CurrentFlux[BinObj_No])
                                        x_flux=x_flux
                                        Plate_X.append(x_flux)
                                        Plate_redshifts.append(CurrentSup.z)
                                        Plate_AND.append(ObjAndM)
                                        Plate_Inv.append(ObjInv)
                                        Plate_Mag.append(CurrentSup.Mag)
                                        Plate_Name.append(CurrentSup.name)
                                        
                                                
                                    
                                else: 
                                    if CurrentSup.Class_p ==1:
                                        Plate_Y.append(0)
                                        x_flux =(CurrentFlux[BinObj_No])
                                        x_flux=x_flux
                                        Plate_X.append(x_flux)
                                        Plate_redshifts.append(CurrentSup.z)
                                        Plate_AND.append(ObjAndM)
                                        Plate_Inv.append(ObjInv)
                                        Plate_Mag.append(CurrentSup.Mag)
                                        Plate_Name.append(CurrentSup.name)
                                        
                                    else:
                                        Plate_Y.append(2)
                                        x_flux =(CurrentFlux[BinObj_No])
                                        x_flux=x_flux
                                        Plate_X.append(x_flux)
                                        Plate_redshifts.append(CurrentSup.z)
                                        Plate_AND.append(ObjAndM)
                                        Plate_Inv.append(ObjInv)
                                        Plate_Mag.append(CurrentSup.Mag)
                                        Plate_Name.append(CurrentSup.name)
                                    
                                        
            
                BinObj_No=BinObj_No+1  
            Sup_obj=Sup_obj+1
        # Only append plates with matched objects
        if plate_match:
            All_Y.append(Plate_Y)
            All_X.append(Plate_X)
            All_Name.append(Plate_Name)
            wav_logs.append(wav)
            All_redshifts.append(Plate_redshifts)
            All_Mag.append(Plate_Mag)
            All_AND.append(Plate_AND)
            All_Inv.append(Plate_Inv)
            All_MJDs.append(currentmjd)
            MatchedPlates.append(pid)
            plate_no=plate_no+1
        else:
            plate_no=plate_no+1
            


    return All_X,All_Y,All_redshifts,All_Mag,All_AND,All_Inv,wav_logs,All_Name,All_MJDs,MatchedPlates


def classification(objectclass, Trainingclass, prediction):
    ## Star, Quasar, Galaxy, BAL
    classi = []
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
    # Checking the classification of each predicted object
    while i< len(Trainingclass):
        currentobject = Trainingclass[i]
        currentpred = prediction[i]
        if objectclass==currentobject:
            if currentpred==0: 
                star=star+1
                starloc.append(i)
            elif currentpred==1: 
                qso=qso+1
                qsoloc.append(i)
            elif currentpred==2: 
                gal=gal+1
                galloc.append(i)
            elif currentpred==3: 
                bal=bal+1
                balloc.append(i)
        
        i=i+1
    
    classi.append(star/(star+qso+gal+bal))
    classi.append(qso/(star+qso+gal+bal))
    classi.append(gal/(star+qso+gal+bal))
    classi.append(bal/(star+qso+gal+bal))
    return classi,starloc,qsoloc,galloc,balloc

def storing(PLATEIDs,supers):
    """Input: Plate IDs and Superset data
       Output: Array of objects with their relevant data
    """
    Full_Data=[]
    matched_objects=0
    print("Number of plates in superset",len(supers))
    for Current_Plate in PLATEIDs:
        platename_data = []
        for Current in supers:
    
            ##Checking if object was surveyed by current plate of interest
            if Current['PLATE'] == Current_Plate:
                    
                
                 ##Storing object data
                Object_ = Object(Current['SDSS_NAME'], Current['RA'], Current['Dec'], 
                         Current['Z_VI'], Current['CLASS_PERSON'],Current['PLATE'] ,
                         Current['MJD'], Current['FIBERID'],Current['PSFMAG'])
                ##Adding object data to array containing data from other objects of the same plate
                platename_data.append(Object_) 
                matched_objects+=1
           
        Full_Data.append(platename_data)
    return Full_Data




