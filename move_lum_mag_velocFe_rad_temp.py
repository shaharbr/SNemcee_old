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
              for Mzams in [9.0]
              for E_final in [0.5, 0.9, 1.3]
              for Ni_mass in [0.02, 0.07, 0.12]
              for Ni_boundary in [2.0, 8.0]
              for K_CSM in [10, 30, 60]
              for R_CSM in [500, 1000, 2000]
              ]


# for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]
# for E_final in [0.5, 0.9, 1.3, 1.7]
# for Ni_mass in [0.02, 0.07, 0.12]
# for Ni_boundary in [2.0, 8.0]
# for K_CSM in [0, 10, 30, 60]
# for R_CSM in [0, 500, 1000, 2000]

file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/lum_observed.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_lum_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_lum_data/' + name + '/lum_observed.dat') for name in list_names]

for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])

file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/vel_photo.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_veloc_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_veloc_data/' + name + '/vel_Fe.dat') for name in list_names]

for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])

file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/magnitudes.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_mag_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_mag_data/' + name + '/magnitudes.dat') for name in list_names]

for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])

dst_dir = [str('/home/sbracha/SNEC/all_temp_rad_data/' + name) for name in list_names]
T_file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/T_eff.dat') for name in list_names]
T_dst_paths = [str('/home/sbracha/SNEC/all_temp_rad_data/' + name + '/T_eff.dat') for name in list_names]
R_file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/rad_photo.dat') for name in list_names]
R_dst_paths = [str('/home/sbracha/SNEC/all_temp_rad_data/' + name + '/rad_photo.dat') for name in list_names]

for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(T_file_paths[i], T_dst_paths[i])
    shutil.copy(R_file_paths[i], R_dst_paths[i])

