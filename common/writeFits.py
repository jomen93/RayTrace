import matplotlib.pyplot as plt
from astropy.io import fits 

class FITS:
    def __init__(self,info,name):
        self.info = info
        self.name = name
        self.hdu = fits.PrimaryHDU(data = self.info)
        
    def Write(self):
        self.hdu.writeto(self.name)
        self.hdu.header['i'] = '1'
    
    def showImage(self):
        plt.imshow(fits.open(self.name)[0].data)
        plt.colorbar()
        plt.show()


        
        
