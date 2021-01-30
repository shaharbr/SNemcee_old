import numpy as np



def distance_luminosity_to_modulus(dL, dL_err):
    dm = 5 * np.log10(dL) - 5
    dm_err = 5 / (np.log(10) * dL) * dL_err
    return dm, dm_err


# dL = 3.43 * 10**6
# dL_err = 0.21 * 10**6
# dm, dm_err = distance_luminosity_to_modulus(dL, dL_err)
# print(dm)
# print(dm_err)


'''
# dL values from:
1) "The distance, supernova rate, and supernova progenitors of NGC 6946" J J Eldridge, Lin Xiao 2019

'''




from astropy import cosmology
from astropy.coordinates import Distance
cosmology.default_cosmology.set(cosmology.WMAP7)
redshift = 0.00167
distmod = Distance(z=redshift).distmod.value
print(distmod)