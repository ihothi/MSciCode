import sys
import numpy as np


b=10
i=0

while i<10:
    
    print()
    sys.stdout.write("\r" +np.str(i*100/10)+"%")
    sys.stdout.flush()
    sys.stdout.flush()
    i=i+1