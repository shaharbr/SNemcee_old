import multiprocessing
import time


def cpu_bound(number):
    return sum(i * i for i in range(number))


def find_sums(numbers):
    with multiprocessing.Pool() as pool:
        pool.map(cpu_bound, numbers)




import asyncio, os
import shutil
import tarfile
import compute_velocity as cv
import datetime

project_dir = os.path.join('/home','sbracha','SNEC','snecmachine_test')
time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def progress_txt(str):
    with open(os.path.join(project_dir, time_now+'_asyncio_w_wind_output.txt'), "a") as text_file:
        text_file.write(str+'\n')


def extract_tar_gz(filename, dst_dir):
    tar = tarfile.open(filename+'.tar.gz', "r:gz")
    tar.extractall(dst_dir)
    tar.close()


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def veloc(model_values):
    [Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore] = model_values
    # only enter if semaphore can be acquired
    async with semaphore:
        if K_CSM == 0 or R_CSM == 0:
            K_CSM = 0
            R_CSM = 0
        name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
               + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
        if E_final < 0.5 and Mzams > 12.0:
            print('skipped impossible SN ' + name)
        elif E_final > 0.5 and Ni_mass < 0.01:
            print('skipped unnecessary SN ' + name)
        model_dir = os.path.join(project_dir, name)
        dat_dir_name = os.path.join(model_dir, name + '_dat')
        xg_dir_name = os.path.join(model_dir, name + '_xg')
        if os.path.exists(dat_dir_name) or os.path.exists(xg_dir_name):
            print('skipped ' + name + ', already exists')
        else:
            print('start ' + name)
            extract_tar_gz(dat_dir_name, model_dir)
            extract_tar_gz(xg_dir_name, model_dir)
            cv.compute_vel(Mzams, dat_dir_name, xg_dir_name, project_dir)
            make_tarfile(dat_dir_name + '.tar.gz', dat_dir_name)
            make_tarfile(xg_dir_name + '.tar.gz', xg_dir_name)
            shutil.rmtree(dat_dir_name)
            shutil.rmtree(xg_dir_name)
            print('end ' + name)



# async def main():
#     semaphore = asyncio.BoundedSemaphore(64)
#     await asyncio.gather(*[veloc(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore)
#                            for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]
#                            for E_final in [0.5, 0.9, 1.3, 1.7]
#                            for Ni_mass in [0.02, 0.07, 0.12]
#                            for Ni_boundary in [2.0, 8.0]
#                            for K_CSM in [0, 10, 30, 60]
#                            for R_CSM in [0, 500, 1000, 2000]
#                            ])

# async def main():
#     semaphore = asyncio.BoundedSemaphore(16)
#     await asyncio.gather(*[veloc(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore)
#                            for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]
#                            for E_final in [0.5, 0.9, 1.3, 1.7]
#                            for Ni_mass in [0.02, 0.07, 0.12]
#                            for Ni_boundary in [2.0, 8.0]
#                            for K_CSM in [0, 10, 30, 60]
#                            for R_CSM in [0, 500, 1000, 2000]
#                            ])




if __name__ == "__main__":
    numbers = [5_000_000 + x for x in range(20)]
    semaphore = multiprocessing.BoundedSemaphore(16)
    models = [[Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore]
    for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]
    for E_final in[0.5, 0.9, 1.3, 1.7]
    for Ni_mass in[0.02, 0.07, 0.12]
    for Ni_boundary in[2.0, 8.0]
    for K_CSM in[0, 10, 30, 60]
    for R_CSM in[0, 500, 1000, 2000]]
    veloc(models)



asyncio.run(main())
