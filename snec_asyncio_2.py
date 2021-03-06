import asyncio, os
import shutil, errno
import time
import re
from add_wind_local_addedR import add_wind


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


async def snec(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore):
    # only enter if semaphore can be acquired
    async with semaphore:
        name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
               + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
        # print('start ' + name)

        dir_name = '/home/sbracha/SNEC/' + name + '/'

        snec_src = '/home/sbracha/SNEC/SNEC_3/'
        copyanything(snec_src, dir_name)

        profile_src = '/home/sbracha/SNEC/sukhbold_profiles/s'+str(int(Mzams))+'/profiles/'
        copyanything(profile_src, dir_name+'/profiles/')
        time.sleep(2)

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

        cmd = './snec > '+name+'.txt'

        proc = await asyncio.create_subprocess_shell(cmd, stdout=None, stderr=None)

        await proc.communicate()
        data_src = dir_name+'Data/'
        shutil.copyfile(dir_name+name+'.txt', data_src+name+'.txt')
        purge(data_src, '.*.xg')
        data_dst = '/home/sbracha/SNEC/all_data/' + name
        copyanything(data_src, data_dst)
        shutil.rmtree(dir_name)

async def main():
    semaphore = asyncio.BoundedSemaphore(16)
    await asyncio.gather(*[snec(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore)
                           for E_final in [0.5, 0.9, 1.3, 1.7]
                           for Mzams in [9.0, 11.0,  13.0, 15.0, 17.0]
                           for Ni_mass in [0.02, 0.12]
                           for Ni_boundary in [2.0, 8.0]
                           for R_CSM in [1500]
                           for K_CSM in [50]
                           ])


asyncio.run(main())
