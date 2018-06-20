

import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits 
from ClassFits import FITS
from DISK import Disk

# Datos de salida de ejemplo
R_final = 50*np.random.rand(40,40)


complete_disk = Disk(R_final)
disk = complete_disk.Ch()
#print disco

name = "example.fits"
image = FITS(disk, name)
image.Write()
image.Pint()

#Info = fits.open(name)
#Info2 = Info[0]
#print Info2.header  ##['i']


#### Falta reescribir el documento y confirmar el header
