import numpy as np
from astropy.io import fits
import os
import desitarget
from desitarget import desi_mask
from desiutil.bitmask import BitMas

slash =  os.path.normpath("/")
#Where will the fules be saved?
desidir=slash+"share"+slash+"splinter"+slash+"hothi"+slash+"Newdesi_sim_data"
#location of DESI data
Platedir = os.path.normpath(slash+"share"+slash+"data1"+slash+"desi_sim"+slash+"spectro"+slash+"redux"+slash+"dc17a2"+slash+"spectra-64")
folder = os.listdir(Platedir)
#Truth File location
Superpath = os.path.normpath(slash+"share"+slash+"data1"+slash+"desi_sim"+slash+"targets"+slash+"truth-lite.fits")

file="folder_Name.txt"
p = open(file, 'w')
for i in folder:
    p.write(i +'\n')

p.close()


#Opening HEALPix data
with open(file) as f:
    folders = f.read().splitlines()

print(np.str(len(folders)))
Folder_all=[]
for i in folders:
    Folder_all.append(i)

groups =[]
FluxB=[]
FluxR=[]
FluxZ=[]
ANDMASKB=[]
ANDMASKR=[]
ANDMASKZ=[]
INVARB=[]
INVARR=[]
INVARZ = []
WAVB = []
WAVR = []
WAVZ = []
ID =[]
print("Opening Files")
for f in Folder_all:
    folder_list = os.listdir(Platedir+slash+f)
    for l in folder_list:
        filefolder_list = os.listdir(Platedir+slash+f+slash+l)
        for n in filefolder_list:
            if 'spectra' in n:
                c=Platedir+slash+f+slash+l+slash+n
                plate_ = fits.open(c,memmap=True)
                wqso =plate_[1].data['DESI_TARGET'][:] & desi_mask.mask('QSO')
                wqso  = wqso >0
                Bin_info_ = plate_[1].data[wqso]
                Flux_B = plate_[3].data[wqso]
                Flux_R = plate_[8].data[wqso]
                Flux_Z = plate_[13].data[wqso]
                AND_B = plate_[5].data[wqso]
                AND_R = plate_[10].data[wqso]
                AND_Z = plate_[15].data[wqso]
                INV_B = plate_[4].data[wqso]
                INV_R = plate_[9].data[wqso]
                INV_Z = plate_[14].data[wqso]
                WAV_B = plate_[2].data
                WAV_R = plate_[7].data
                WAV_Z = plate_[12].data
                FluxB.append(Flux_B)
                FluxR.append(Flux_R)
                FluxZ.append(Flux_Z)
                groups.append(n)
                ANDMASKB.append(AND_B)
                INVARB.append(INV_B)
                ANDMASKR.append(AND_R)
                INVARR.append(INV_R)
                ANDMASKZ.append(AND_Z)
                INVARZ.append(INV_Z)
                WAVB.append(WAV_B)
                WAVR.append(WAV_R)
                WAVZ.append(WAV_Z)
                ID.append(Bin_info_['TARGETID'])
  

#Opening Truth File      
Super = fits.open(Superpath,memmap=True)
TargetID = Super[1].data['TARGETID']
Z = Super[1].data['TRUEZ']
SPECTYPE = Super[1].data['TRUESPECTYPE']

match_name = []
match_FluxB=[]
match_FluxR=[]
match_FluxZ=[]
match_ANDMASKB=[]
match_ANDMASKR=[]
match_ANDMASKZ=[]
match_INVARB=[]
match_INVARR=[]
match_INVARZ = []
match_WAVB =[]
match_WAVR = []
match_WAVZ = []
match_ID =[]
match_Z = []
match_Class =[]




total_match=0
g=0


#Matching HEALPix objects to truth file

while g <len(groups):
    matched = False 
    m_FluxB=[]
    m_FluxR=[]
    m_FluxZ=[]
    m_ANDMASKB=[]
    m_ANDMASKR=[]
    m_ANDMASKZ=[]
    m_INVARB=[]
    m_INVARR=[]
    m_INVARZ = []
    m_WAVB =[]
    m_WAVR = []
    m_WAVZ = []
    m_ID =[]
    m_Z =[]
    m_Class =[]

    plate_name = groups[g]
    plate_FluxB=FluxB[g]
    plate_FluxR=FluxR[g]
    plate_FluxZ=FluxZ[g]
    plate_ANDMASKB=ANDMASKB[g]
    plate_ANDMASKR=ANDMASKR[g]
    plate_ANDMASKZ=ANDMASKZ[g]
    plate_INVARB=INVARB[g]
    
    plate_INVARR=INVARR[g]
    plate_INVARZ = INVARZ[g]
    plate_WAVB = WAVB[g]
    plate_WAVR = WAVR[g]
    plate_WAVZ = WAVZ[g]
    plate_ID =ID[g]
    n=0
    while n<len(plate_FluxB):
        current_FluxB=plate_FluxB[n]
        current_FluxR=plate_FluxR[n]
        current_FluxZ=plate_FluxZ[n]
        current_ANDMASKB=plate_ANDMASKB[n]
        current_ANDMASKR=plate_ANDMASKR[n]
        current_ANDMASKZ=plate_ANDMASKZ[n]
        current_INVARB=plate_INVARB[n]
        current_INVARR=plate_INVARR[n]
        current_INVARZ = plate_INVARZ[n]
        current_ID =plate_ID[n]
        s=0
        while s< len(TargetID):
            super_id = TargetID[s]
            if current_ID==super_id:
                matched=True
                m_FluxB.append(current_FluxB)
                m_FluxR.append(current_FluxR)
               
                
                m_FluxZ.append(current_FluxZ)
                m_ANDMASKB.append(current_ANDMASKB)
                m_ANDMASKR.append(current_ANDMASKR)
                m_ANDMASKZ.append(current_ANDMASKZ)
                m_INVARB.append(current_INVARB)
                m_INVARR.append(current_INVARR)    
                m_INVARZ.append(current_INVARZ)
                m_ID.append(current_ID)
                m_Z.append(Z[s])
                if SPECTYPE[s] =='GALAXY':
                    m_Class.append(2)
                if SPECTYPE[s] =='STAR':
                    m_Class.append(0)
                if SPECTYPE[s] =='QSO':
                    if Z[s]<2.1:
                        
                        m_Class.append(1)
                    else:
                        m_Class.append(3)
            s+=1
        n+=1
        
    if matched:
        match_name.append(plate_name)
        match_FluxB.append(m_FluxB)
        match_FluxR.append(m_FluxR)
        match_FluxZ.append(m_FluxZ)
        match_ANDMASKB.append(m_ANDMASKB)
        match_ANDMASKR.append(m_ANDMASKR)
        match_ANDMASKZ.append(m_ANDMASKZ)
    
        match_INVARB.append(m_INVARB)
        match_INVARR.append(m_INVARR)
        match_INVARZ.append(m_INVARZ)
        match_WAVB.append(plate_WAVB)
        match_WAVR.append(plate_WAVR)
        match_WAVZ.append(plate_WAVZ)
        match_ID.append(m_ID)
        match_Z.append(m_Z)
        match_Class.append(m_Class)        
    g+=1
    








i=0
#Saving files
First_Run=True
while i <len(match_FluxB):
    C_Fiber = match_name[i]
    a1=match_FluxB[i]
    a2=match_ANDMASKB[i] 
    a3=match_INVARB[i]
    print(a3)
    a4=match_WAVB[i]
    a5=match_FluxR[i]
    a6=match_ANDMASKR[i] 
    a7=match_INVARR[i]
    a8=match_WAVR[i]
    a9=match_FluxZ[i]
    a10=match_ANDMASKZ[i] 
    a11=match_INVARZ[i]
    a12=match_WAVZ[i]
    a13=match_Class[i]
    a14=match_Z[i]
    a15=match_ID[i]
    col1 = fits.Column(name='BFLUX', format='PD()', array=np.array(a1,dtype=np.object))
    col2 = fits.Column(name='BAND', format='PI()', array=np.array(a2,dtype=np.object))
    col3 = fits.Column(name='BINV', format='PD()', array=np.array(a3,dtype=np.object))
    col4 = fits.Column(name='BWAV', format='D', array=np.array(a4))
    col5 = fits.Column(name='RFLUX', format='PD()', array=np.array(a5,dtype=np.object))
    col6 = fits.Column(name='RAND', format='PI()', array=np.array(a6,dtype=np.object))
    col7 = fits.Column(name='RINV', format='PD()', array=np.array(a7,dtype=np.object))
    col8 = fits.Column(name='RWAV', format='D', array=np.array(a8))
    col9 = fits.Column(name='ZFLUX', format='PD()', array=np.array(a9,dtype=np.object))
    col10 = fits.Column(name='ZAND', format='PI()', array=np.array(a10,dtype=np.object))
    col11 = fits.Column(name='ZINV', format='PD()', array=np.array(a11,dtype=np.object))
    col12 = fits.Column(name='ZWAV', format='D', array=np.array(a12))
    col13 = fits.Column(name='SPECTYPE', format='I', array=np.array(a13))
    col14 = fits.Column(name='TRUEZ', format='D', array=np.array(a14))
    col15 = fits.Column(name='TARGETID', format='I', array=np.array(a15))
    cols = fits.ColDefs([col1, col2,col3, col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr['FIBER'] = C_Fiber
    prihdr['OBJNO'] = len(a5)
    prihdu = fits.PrimaryHDU(header=prihdr)
    file_dir = desidir+slash+np.str(C_Fiber)+slash
    os.mkdir(file_dir)
    file_name = file_dir+"New"+np.str(C_Fiber)
    print(file_name)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(file_name)
    i=i+1
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    