import numpy as np
from astropy.io import fits

def Reading(list):
    """ Input: list - name of the .txt file that contains list of .fits files 
        Output: F - a list containing names of .fits files"""
    F = open(list, 'w') 
    return F


class Object:
    
    def __init__():[]
   
    def __init__(self, SDSS_NAME, R, D, Z_VI, CLASS_P, p, mjd, fid):
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
        