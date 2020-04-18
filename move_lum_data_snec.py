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
              for Mzams in [16.0]
              for Ni_mass in [0.10]
              for E_final in [1.0, 1.4, 1.8, 2.2]
              for Ni_boundary in [3.0]
              for R_CSM in [500]
              for K_CSM in [0.001, 20, 40, 60, 80]
              ]

file_paths = [str('/home/sbracha/SNEC/all_data/' + name + '/lum_observed.dat') for name in list_names]
dst_dir = [str('/home/sbracha/SNEC/all_lum_data/' + name) for name in list_names]
dst_paths = [str('/home/sbracha/SNEC/all_lum_data/' + name + '/lum_observed.dat') for name in list_names]

list_names = ['dir1', 'dir2']

# file_paths = [str('/home/sbracha/test/' + name + '/file1.txt') for name in list_names]
# dst_dir = [str('/home/sbracha/test/all/' + name) for name in list_names]
# dst_paths = [str('/home/sbracha/test/all/' + name + '/file1.txt') for name in list_names]


for i in range(len(list_names)):
    os.mkdir(dst_dir[i])
    shutil.copy(file_paths[i], dst_paths[i])

