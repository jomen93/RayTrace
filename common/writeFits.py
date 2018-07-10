import matplotlib.pyplot as plt
from astropy.io import fits 

class FITS:
    def __init__(self,info,name, hdrData):
        self.info = info
        self.name = name
        self.hdrData = hdrData
        self.hdu = fits.PrimaryHDU(data = self.info)

        
    def Write(self):
        self.hdu.writeto(self.name)
        self.hdu.header['i'] = '1'

        # Include data of the image in the Header of the .fits file
        for j in self.hdrData.keys():
            fits.setval(self.name, j, value=self.hdrData[j])

    def showImage(self):
        plt.imshow(fits.open(self.name)[0].data)
        plt.colorbar()
        plt.show()


        
        
