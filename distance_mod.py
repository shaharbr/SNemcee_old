from astropy import cosmology
from astropy.coordinates import Distance
cosmology.default_cosmology.set(cosmology.WMAP7)
redshift = 0.025
distmod = Distance(z=redshift).distmod.value
print(distmod)