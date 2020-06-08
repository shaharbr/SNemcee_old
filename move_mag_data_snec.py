import os
import shutil, errno


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


list_names = [str('M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)\
                   + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM))
              for Mzams in [13.0, 16.0, 19.0, 21.0]
              for Ni_mass in [0.10, 0.13, 0.16, 0.19]
              for E_final in [1.5, 1.8, 2.4]
              for Ni_boundary in [3.0]
              for R_CSM in [600, 1800, 2400, 3000]
              for K_CSM in [0.001, 30, 90]
              ]

file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/lum_observed.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_mag_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_mag_data/' + name + '/lum_observed.dat') for name in list_names]


for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])

