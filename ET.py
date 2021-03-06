import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from sklearn.metrics import auc


def calculate_ET(time_vec, Lum_vec, Ni_vec):
    # integ_Lum = auc(time_vec, time_vec * Lum_vec)
    # integ_Ni = auc(time_vec, time_vec * Ni_vec)
    # ET = auc(time_vec, time_vec * (Lum_vec - Ni_vec))
    function = [time_vec[i] * (Lum_vec[i] - Ni_vec[i]) for i in range(len(Lum_vec))]
    ET = [auc(time_vec[:i], function[:i]) for i in np.arange(2, len(Lum_vec)+1)]
    return ET


blackbody_data = pd.read_csv(r'results\blackbody_results_18hmx.csv')

Lum_vec = list(blackbody_data['Lum']) #  erg/s
Lum_vec = [x * 86400 * 10**7 for x in Lum_vec] #  erg/s
# print(Lum_vec)


time_vec = list(blackbody_data['t_from_discovery']) #  erg/s
print(time_vec)
# print(time_vec)

len = len(time_vec)

Ni_mass = pd.read_csv(r'results\Ni_mass_SN2018hmx.csv')

Ni_vec = list(Ni_mass['Ni_mass']) #  erg/s
Ni_vec = [x * 86400 * 10**43 for x in Ni_vec] #  erg/s
# print(Ni_vec)


ET_vec = calculate_ET(time_vec, Lum_vec, Ni_vec)
# print(ET_vec)

plt.figure()
plt.plot(time_vec, [0]+ET_vec, label='2018hmx', marker='.')
plt.ylabel('ET (erg/days)')
plt.xlabel('rest days from discovery')
plt.title('ET = ∫ L tdt-∫ QNi tdt [erg d]')
plt.xlim(left=0)
plt.legend()


plt.figure()
plt.plot(time_vec, [0]+[ET_vec[i] * 86400 for i in range(len-1)], label='2018hmx', marker='.')
plt.ylabel('ET (erg/s)')
plt.xlabel('rest days from discovery')
plt.title('ET = ∫ L tdt-∫ QNi tdt [erg s]')
plt.xlim(left=0)
plt.legend()
plt.ylim(bottom=0)

Fe_v = 5600 * (10**5) # cm/s


final_ET =  ET_vec[-1]  * 86400# erg/s -> g*cm^2/s^3
print(final_ET)

# ET/v50
# (g * cm^2 / s^3) / (cm/s) = g * cm
ET_div_v50 = final_ET / Fe_v
print(ET_div_v50)

# TODO save as csv

plt.show()