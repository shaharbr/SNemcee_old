import pandas as pd
import bolometric

from matplotlib import pyplot as plt

blackbody_data = pd.read_csv('blackbody_results.csv')

Lum_vec = list(blackbody_data['Lum']) #  erg/s
Lum_vec = [x * 86400 * 10**7 for x in Lum_vec] #  erg/s
# print(Lum_vec)


time_vec = list(blackbody_data['t_from_discovery']) #  erg/s
# print(time_vec)

len = len(time_vec)

Ni_mass = pd.read_csv('Ni_mass_SN2018hmx.csv')

Ni_vec = list(Ni_mass['Ni_mass']) #  erg/s
Ni_vec = [x * 86400 * 10**43 for x in Ni_vec] #  erg/s
# print(Ni_vec)


ET_vec = bolometric.calculate_ET(time_vec, Lum_vec, Ni_vec)
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

print([ET_vec[i] * 86400 for i in range(len-2)])

plt.show()