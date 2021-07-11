import multiprocessing
# import asyncio
import os
import shutil, errno
import tarfile
import compute_velocity as cv
import datetime
import time

project_dir = os.path.join('/home','sbracha','SNEC')
outputs_dir = os.path.join(project_dir, 'all_snec_outputs')
time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise


def progress_txt(str):
    with open(os.path.join(outputs_dir, time_now+'_asyncio_w_wind_output.txt'), "a") as text_file:
        text_file.write(str+'\n')


def extract_tar_gz(filename, dst_dir):
    tar = tarfile.open(filename+'.tar.gz', "r:gz")
    tar.extractall(dst_dir)
    tar.close()


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def veloc(Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM, semaphore):
    # only enter if semaphore can be acquired
    if K_CSM == 0 or R_CSM == 0:
        K_CSM = 0
        R_CSM = 0
    name = 'M' + str(Mzams) + '_Ni' + str(Ni_mass) + '_E' + str(E_final) \
           + '_Mix' + str(Ni_boundary) + '_R' + str(R_CSM) + '_K' + str(K_CSM)
    if E_final < 0.5 and Mzams > 12.0:
        print('skipped impossible SN ' + name)
    elif E_final > 0.5 and Ni_mass < 0.01:
        print('skipped unnecessary SN ' + name)
    model_dir = os.path.join(outputs_dir, name)
    dat_dir_name = os.path.join(model_dir, name + '_dat')
    xg_dir_name = os.path.join(model_dir, name + '_xg')
    print('start ' + name)
    extract_tar_gz(dat_dir_name, model_dir)
    time.sleep(30)
    veloc_file_path = os.path.join(dat_dir_name, 'vel_Fe.dat')
    if os.path.exists(veloc_file_path):
        print('skipped ' + name + ', vel_Fe.dat already exists')
        make_tarfile(dat_dir_name + '.tar.gz', dat_dir_name)
        shutil.rmtree(dat_dir_name)
        sema.release()
    else:
        extract_tar_gz(xg_dir_name, model_dir)
        time.sleep(30)
        cv.compute_vel(Mzams, dat_dir_name, xg_dir_name, project_dir)
        make_tarfile(dat_dir_name + '.tar.gz', dat_dir_name)
        make_tarfile(xg_dir_name + '.tar.gz', xg_dir_name)
        dat_dir_dst = os.path.join(project_dir, 'all_data', name)
        os.rename(dat_dir_name, dat_dir_dst)
        shutil.rmtree(xg_dir_name)
        print('end ' + name)
        sema.release()


if __name__ == '__main__':
    concurrency = 20
    sema = multiprocessing.Semaphore(concurrency)
    all_processes = []

    models = [[Mzams, Ni_mass, E_final, Ni_boundary, R_CSM, K_CSM]
              for Mzams in [9.0, 11.0, 13.0, 15.0, 17.0]
              for E_final in [0.5, 0.9, 1.3, 1.7]
              for Ni_mass in [0.02, 0.07, 0.12]
              for Ni_boundary in [1.41]
              for K_CSM in [500, 1000, 2000]
              for R_CSM in [10, 30, 60]]

    for model in models:
        # once 20 processes are running, the following `acquire` call
        # will block the main process since `sema` has been reduced
        # to 0. This loop will continue only after one or more
        # previously created processes complete.
        sema.acquire()
        p = multiprocessing.Process(target=veloc, args=(*model, sema))
        all_processes.append(p)
        p.start()

    # inside main process, wait for all processes to finish
    for p in all_processes:
        p.join()
