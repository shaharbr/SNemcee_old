import asyncio, os
import shutil, errno
import time
import re

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


async def snec(Mzams, Ni_mass, E_final, Ni_boundary, semaphore):
    # only enter if semaphore can be acquired
    async with semaphore:
        name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final)\
                   + '_Mix' + str(Ni_boundary) + '_noCSM'
        # print('start ' + name)

        dir_name = '/home/sbracha/' + name + '/'

        snec_src = '/home/sbracha/SNEC_template/'
        copyanything(snec_src, dir_name)

        profile_src = '/home/sbracha/sukhbold_profiles/s'+str(Mzams)+'/profiles/'
        copyanything(profile_src, dir_name+'/profiles/')
        time.sleep(2)

        os.chdir(os.path.abspath(dir_name))

        # Read in the file
        with open('parameters', 'r') as file:
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('profiles/xxx.', 'profiles/s'+str(Mzams)+'.')
        filedata = filedata.replace('Ni_mass = xxx', 'Ni_mass = '+str(Ni_mass))
        filedata = filedata.replace('final_energy        = xxxd51', 'final_energy        = ' + str(E_final)+'d51')
        filedata = filedata.replace('Ni_boundary_mass = xxxd0', 'Ni_boundary_mass = ' + str(Ni_boundary) + 'd0')
        # Write the file out again
        with open('parameters', 'w') as file:
            file.write(filedata)

        cmd = './snec > '+name+'.txt'

        proc = await asyncio.create_subprocess_shell(cmd, stdout=None, stderr=None)

        await proc.communicate()

        data_src = dir_name+'/Data/'
        purge(data_src, '.*.xg')
        data_dst = '/home/sbracha/all_data/' + name
        copyanything(data_src, data_dst)
        shutil.rmtree(dir_name)

async def main():
    semaphore = asyncio.BoundedSemaphore(1)
    await asyncio.gather(*[snec(Mzams, Ni_mass, E_final, Ni_boundary, semaphore)
                           for Mzams in [13.0, 15.0]
                           for Ni_mass in [0.02, 0.12]
                           for E_final in [0.9, 1.3]
                           for Ni_boundary in [5.0, 8.0]
                         ])


asyncio.run(main())
