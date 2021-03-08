import glob
import numpy
import re

#import scidata.monodataset
#import scidata.utils

import monodataset
import utils

def parse(liter, column=None):
    """
    Imports the data from a given iterator and returns a
    scidata.monodataset.dataset object

    Example data = xgraph.parse(open("data.xg","r"))
    """
#    dataset = scidata.monodataset.dataset()
    dataset = monodataset.dataset()

    if column is None:
        column = 1
    else:
        column -= 1

    dataset.nframes  = 0
    dataset.time     = []
    dataset.metadata = []
    dataset.data_x   = []
    dataset.data_y   = []
    dataset.ranges   = []

    si = 0
    ei = 0
    for l in liter:
        header = re.match(
                r"\s*[#\"]*\s*(\w*)\s*=\s*(-?\d*\.?\d*[eE]?[-+]?\d*).*", l)

        if header is not None:
            if header.group(1) == 'Time':
                if si != ei:
                    dataset.ranges.append((si,ei))
                si = ei
                dataset.time.append(float(header.group(2)))
                dataset.nframes += 1
                dataset.metadata.append({})
            else:
                dataset.metadata[-1][header.group(1)] = header.group(2)
        elif l[0] not in ('#', '"'):
            dataline = l.split()
            if len(dataline) > column:
                dataset.data_x.append(float(dataline[0]))
                dataset.data_y.append(float(dataline[column]))
                ei += 1
    dataset.ranges.append((si,ei))

    dataset.data_x = numpy.array(dataset.data_x)
    dataset.data_y = numpy.array(dataset.data_y)

    return dataset

def parsefile(filename, column=None):
    """
    Import data from a given file
    """
#    if not scidata.utils.extension(filename) in ["xg", "yg"]:
#        raise scidata.utils.FileTypeError(filename)
        
    if not utils.extension(filename) in ["xg", "yg"]:
        raise utils.FileTypeError(filename)

    return parse(open(filename,"r"), column)

def parsefiles(filelist, column=None):
    """
    Loads all the files in the list

    This function returns a dictionary {filename: dataset}, where
    * filename : is the file name
    * dataset  : is the scidata.monodataset.dataset object
    """
    return dict([(f, parsefile(f, column)) for f in filelist])

def loaddir(directory, column=None):
    """
    Loads all the xgraph data files in the given directory into memory

    This function returns a dictionary { varname: dataset }, where
    * varname : is the variable name (taken from the filename)
    * dataset : is the monodataset.dataset object
    """
    out = []

    for f in glob.glob(directory + "/*.xg"):
        out.append((scidata.utils.basename(f), parsefile(f, column)))
    return dict(out)
