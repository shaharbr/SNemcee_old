import asyncio, os
# import numpy as np
import shutil, errno
from add_wind_local import add_wind


# for Mzams in np.arange(11, 21, 0.5):
#     for Ni_mass in np.arange(0.05, 0.21, 0.02):
#         for E_final in np.arange(0.5, 3, 0.2):
#             for Ni_boundary in np.arange(3, Mzams, 1):
#                 for R_CSM in np.arange(500, 5000, 500):
#                     for K_CSM in np.arange(0.01, 100, 10):
# for E_final in [1, 2]:
            #     for Ni_boundary in [3]:
            #         for R_CSM in [1000, 2000]:
            #             for K_CSM in [10, 50]:
    # copy SNEC dir
    # dir_name = '/home/sbracha/SNEC/'+'M'+str(Mzams)+'_Ni'+str(Ni_mass)+'_E'+str(E_final)+'_Mix'+str(Ni_boundary)+'_R'+str(R_CSM)+'_K'+str(K_CSM)


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


async def snec(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore):
    # only enter if semaphore can be acquired
    async with semaphore:
        name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)\
                   + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
        print('start ' + name)

        dir_name = '/home/sbracha/SNEC/' + name + '/'

        snec_src = '/home/sbracha/SNEC/SNEC/'
        copyanything(snec_src, dir_name)

        profile_src = '/home/sbracha/SNEC/sukhbold_profiles/s'+str(int(Mzams))+'/profiles/'
        copyanything(profile_src, dir_name+'/profiles/')

        os.chdir(os.path.abspath(dir_name))

        add_wind(K_CSM, R_CSM, Mzams, dir_name)

        # Read in the file
        with open('parameters', 'r') as file:
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('profiles/xxx.', 'profiles/' + 's{}_K{}_R{}.'.format(str(Mzams), str(K_CSM), str(R_CSM)))
        filedata = filedata.replace('Ni_mass = xxx', 'Ni_mass = '+str(Ni_mass))
        filedata = filedata.replace('final_energy        = xxxd51', 'final_energy        = ' + str(E_final)+'d51')
        filedata = filedata.replace('Ni_boundary_mass = xxxd0', 'Ni_boundary_mass = ' + str(Ni_boundary) + 'd0')
        # Write the file out again
        with open('parameters', 'w') as file:
            file.write(filedata)

        cmd = './snec'
        proc = await asyncio.create_subprocess_shell(
            cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE)
        stdout, stderr = await proc.communicate()

        data_src = dir_name+'/Data/'
        data_dst = '/home/sbracha/SNEC/all_data/' + name
        copyanything(data_src, data_dst)
        # os.rmdir(dir_name)

        print(f'[{cmd!r} exited with {proc.returncode}]')
        if stdout:
            print(f'[stdout]\n{stdout.decode()}')
        if stderr:
            print(f'[stderr]\n{stderr.decode()}')
        print('end ' + name)


async def main():
    semaphore = asyncio.BoundedSemaphore(12)
    await asyncio.gather(*[snec(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore)
                           for Mzams in [10.0, 13.0, 16.0, 19.0]
                           for Ni_mass in [0.10, 0.13, 0.16]
                           for E_final in [1.0, 1.4, 1.8, 2.2]
                           for Ni_boundary in [3.0]
                           for R_CSM in [500, 1000, 1500, 2000, 2500]
                           for K_CSM in [0.001, 20, 40, 60, 80]
                         ])



# for Mzams in np.arange(11, 21, 0.5):
#     for Ni_mass in np.arange(0.05, 0.21, 0.02):
#         for E_final in np.arange(0.5, 3, 0.2):
#             for Ni_boundary in np.arange(3, Mzams, 1):
#                 for R_CSM in np.arange(500, 5000, 500):
#                     for K_CSM in np.arange(0.01, 100, 10):


asyncio.run(main())
