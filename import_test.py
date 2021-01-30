import os
import shutil

eukdir = os.path.join('/davidb', 'shaharbr','eukProteins_2.5k')

filetypes = {'proteins': ['', '.contigs.min2.5k.proteins.faa', '_proteins.faa'],
             'gff': ['_gff', '.contigs.min2.5k.gff', '_proteins.gff'],
             'dna': ['_dna', '.contigs.min2.5k.fa', '_contigs.fna']}

for source in ['Genomes', 'Metagenomes']:
    source_dir = os.path.join(eukdir, source)
    file_list = os.listdir(source_dir)
    file_list = [f.replace('.contigs.min2.5k.proteins.faa', '')
                        for f in file_list if "all_euk" not in f]
    for filename in file_list:
        dst_dir_path = os.path.join(eukdir, 'all_samples', filename)
        os.mkdir(dst_dir_path)
        for filetype in ['proteins', 'gff', 'dna']:
            datafile_path = os.path.join(eukdir,
                                         source+filetypes[filetype][0],
                                         filename+filetypes[filetype][1])
            distfile_path = os.path.join(dst_dir_path,
                                         filename + filetypes[filetype][2])

            print(filename)
            print(distfile_path)
            print(datafile_path)
            shutil.copy(datafile_path,
                        distfile_path, follow_symlinks=False)
