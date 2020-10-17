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
			for Mzams in [9.0, 18.0]
                        for Ni_mass in [0.02, 0.07, 0.12, 0.17]
                        for E_final in [1.2, 1.7, 2.2, 2.7]
                        for Ni_boundary in [1.0, 3.0, 6.0]
                        for R_CSM in [600, 1400, 2200, 3000]
                        for K_CSM in [0.001, 50, 100, 150]

              ]

file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/lum_observed.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_lum_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_lum_data/' + name + '/lum_observed.dat') for name in list_names]

for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])

file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/vel_photo.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_veloc_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_veloc_data/' + name + '/vel_photo.dat') for name in list_names]


for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])


file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/magnitudes.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_mag_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_mag_data/' + name + '/magnitudes.dat') for name in list_names]


for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])
