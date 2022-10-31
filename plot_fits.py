import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.io import fits
plt.style.use(astropy_mpl_style)
image_data = fits.getdata('/home/Mathieu/Documents/TRAPPIST/reduced_data/CK17K020/20220921TS/images/TRAP.2022-09-22T01:27:24.fits', ext=0)
plt.figure()
plt.imshow(image_data, vmin=0, vmax=60000)
plt.colorbar()
plt.show()
