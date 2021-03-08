import bisect
import copy
import h5py
import math
import numpy
import os
import tempfile
#import scidata.utils
import utils

class frame:
    """
    A single timeframe from a datafile

    Members
    * index    : frame index
    * time     : frame time
    * metadata : dictionary
    * data_x   : x-data array
    * data_y   : y-data array
    """

    def __init__(self):
        """
        Initialize
        """
        self.index = 0
        self.time = 0.0
        self.metadata = {}
        self.data_x = numpy.array([])
        self.data_y = numpy.array([])

    def __call__(self, x):
        """
        Perform a linear interpolation of the data in a point
        """
        return numpy.interp(x, self.data_x, self.data_y)

    def __cmp__(self, other):
        """
        Compare with another dataset
        """
        return cmp(self.time, other.time)

    def __hash__(self):
        return hash(self.time)

    def __lt__(self, other):
        return self.__cmp__(other) < 0
    def __gt__(self, other):
        return self.__cmp__(other) > 0
    def __eq__(self, other):
        return self.__cmp__(other) == 0
    def __le__(self, other):
        return self.__cmp__(other) <= 0
    def __ge__(self, other):
        return self.__cmp__(other) >= 0
    def __ne__(self, other):
        return self.__cmp__(other) != 0

    def __abs__(self):
        out = copy.deepcopy(self)
        out.data_y = abs(out.data_y)
        return out
    def __add__(self, other):
        out = copy.deepcopy(self)
        if isinstance(other, frame):
            out.data_y = self.data_y + other(self.data_x)
        else:
            out.data_y = self.data_y + other
        return out
    __radd__ = __add__
    def __sub__(self, other):
        out = copy.deepcopy(self)
        if isinstance(other, frame):
            out.data_y = self.data_y - other(self.data_x)
        else:
            out.data_y = self.data_y - other
        return out
    def __rsub__(self, other):
        out = copy.deepcopy(self)
        out.data_y = other - self.data_y
        return out
    def __neg__(self):
        out = copy.deepcopy(self)
        out.data_y = - out.data_y
        return out
    def __mul__(self, other):
        out = copy.deepcopy(self)
        if isinstance(other, frame):
            out.data_y = self.data_y * other(self.data_x)
        else:
            out.data_y = self.data_y * other
        return out
    __rmul__ = __mul__
    def __div__(self, other):
        out = copy.deepcopy(self)
        if isinstance(other, frame):
            out.data_y = self.data_y / other(self.data_x)
        else:
            out.data_y = self.data_y / other
        return out
    def __rdiv__(self, other):
        out = copy.deepcopy(self)
        out.data_y = other / self.data_y
        return out
    def __pow__(self, other):
        out = copy.deepcopy(self)
        out.data_y = self.data_y**other
        return out

    def array(self):
        """
        Convert to a 2D array
        """
        return numpy.array((self.data_x, self.data_y))

    def merge(self, other):
        """
        Merge the data from another frame
        """
        sdata_x = self.data_x
        sdata_y = self.data_y
        odata_x = other.data_x
        odata_y = other.data_y

        slist = [(sdata_x[i], sdata_y[i]) for i in range(len(sdata_x))]
        olist = [(odata_x[i], odata_y[i]) for i in range(len(odata_x))]

        tlist = slist + olist
        tlist = sorted(tlist)

        self.data_x = numpy.array([x[0] for x in tlist])
        self.data_y = numpy.array([x[1] for x in tlist])

    def format(self):
        """
        Convert to string in gnuplot format
        """
        s = ['#Time = ' + str(self.time) + '\n']
        for k, m in self.metadata.iteritems():
            s.append('#{0} = {1}\n'.format(k, m))
        return s + ["%.19g %.19g\n" % (self.data_x[i], self.data_y[i]) \
             for i in xrange(self.data_x.size)]

    def plot(self):
        """
        Visualize using pyplot
        """
        import matplotlib.pyplot
        return matplotlib.pyplot.plot(self.data_x, self.data_y)

    def purge_nans(self):
        """
        Remove all the NaNs
        """
        self.data_y = numpy.ma.masked_invalid(self.data_y)
        self.data_x = numpy.ma.masked_array(self.data_x,
                mask=self.data_y.mask)

    def sort(self):
        """
        Sorts the frame
        """
        L = [(self.data_x[i], self.data_y[i]) for i in range(len(self.data_x))]
        L = numpy.array(sorted(set(L)))

        self.data_x = L[:,0]
        self.data_y = L[:,1]

    def write(self,filename):
        """
        Export the data in gnuplot format
        """
        file = open(filename,"w")
        file.writelines(self.format())
        file.close()

class dataset:
    """
    1D database class

    Members:
    * nframes  : frame indices
    * time     : frame times
    * metadata : list of metadata dictionaries
    * data_x   : x-data array for all the frames
    * data_y   : y-data array for all the frames
    * ranges   : indices for the slicing of the data arrays
    """

    def __init__(self):
        self.nframes = 0
        self.time = []
        self.metadata = []
        self.data_x = numpy.array([])
        self.data_y = numpy.array([])
        self.ranges = []

    def __iter__(self):
        return iterator(self)

    def __add__(self, other):
        L = []
        for frame in self:
            if isinstance(other, dataset):
                L.append(frame + other.find_frame(frame.time))
            else:
                L.append(frame + other)
        out = dataset()
        out.import_framelist(L)
        return out
    __radd__ = __add__
    def __sub__(self, other):
        L = []
        for frame in self:
            if isinstance(other, dataset):
                L.append(frame - other.find_frame(frame.time))
            else:
                L.append(frame - other)
        out = dataset()
        out.import_framelist(L)
        return out
    def __rsub__(self, other):
        L = []
        for frame in self:
            L.append(other - frame)
        out = dataset()
        out.import_framelist(L)
        return out
    def __neg__(self):
        L = []
        for frame in self:
            L.append(-frame)
        out = dataset()
        out.import_framelist(L)
        return out
    def __mul__(self, other):
        L = []
        for frame in self:
            if isinstance(other, dataset):
                L.append(frame * other.find_frame(frame.time))
            else:
                L.append(frame * other)
        out = dataset()
        out.import_framelist(L)
        return out
    __rmul__ = __mul__
    def __div__(self, other):
        L = []
        for frame in self:
            if isinstance(other, dataset):
                L.append(frame / other.find_frame(frame.time))
            else:
                L.append(frame / other)
        out = dataset()
        out.import_framelist(L)
        return out
    def __rdiv__(self, other):
        L = []
        for frame in self:
            L.append(other / frame)
        out = dataset()
        out.import_framelist(L)
        return out
    def __pow__(self, other):
        L = []
        for frame in self:
            L.append(frame**other)
        out = dataset()
        out.import_framelist(L)
        return out

    def array(self):
        """
        Convert to a 3D array

        The format will be a_{ijk}
        * i : frame index
        * j : x-data index
        * k : y-data index
        """
        L = []
        for frame in self:
            L.append(frame.array())
        return numpy.array(L)

    def find_frame(self, time):
        """
        Finds the closest frame to the given time
        """
        idx_a = max(bisect.bisect_right(self.time, time)-1, 0)
        if idx_a == self.nframes - 1:
            return self.frame(self.nframes - 1)
        else:
            idx_b = idx_a + 1
            if abs(self.time[idx_a] - time) < abs(self.time[idx_b] - time):
                return self.frame(idx_a)
            else:
                return self.frame(idx_b)

    def frame(self, number):
        """
        Get the wanted frame
        """
        assert number >= 0 and number < self.nframes

        out = frame()

        if self.nframes == 1:
            out.index = 0
            out.time = self.time[0]
            out.data_x = self.data_x
            out.data_y = self.data_y
            out.metadata = self.metadata[0]
        else:
            out.index = number
            out.time = self.time[number]
            out.data_x = self.data_x[self.ranges[number][0]:
                    self.ranges[number][1]]
            out.data_y = self.data_y[self.ranges[number][0]:
                    self.ranges[number][1]]
            out.metadata = self.metadata[number]

        return out

    def format(self):
        """
        Outputs the data in xgraph format

        This will return a string representation of the data
        """
        out = []
        for f in self:
            out.append("\n\n")
            out += f.format()
        return out

    def import_array(self,array):
        """
        Import the data from a full numpy.array

        The array must be a 2D array with the format [n,t,x,y], where
        * n : 0,1,2,3, ... is the frame number
        * t : is the time (real)
        * x : is the x-data array (numpy array)
        * y : is the y-data array (numpy array)
        """
        si = 0
        ei = 0

        self.nframes = int(max(array[:,0]))+1
        self.time = []
        self.data_x = []
        self.data_y = []
        self.metadata = []
        self.ranges = []

        for k in xrange(self.nframes):
            slice = array[array[:,0] == k, 1:4]

            self.time.append(slice[0, 0])

            self.metadata.append({})

            self.data_x += list(slice[:, 1])
            self.data_y += list(slice[:, 2])

            ei += slice[:, 0].size
            self.ranges.append((si,ei))
            si = ei
        self.data_x = numpy.array(self.data_x)
        self.data_y = numpy.array(self.data_y)

    def import_framelist(self,framelist):
        """
        Import the data from a list of frames
        """
        si = 0
        ei = 0

        self.nframes = len(framelist)
        self.time = []
        self.metadata = []
        self.data_x = []
        self.data_y = []
        self.ranges = []

        k = 0
        for f in framelist:
            self.time.append(f.time)

            self.metadata.append(f.metadata)
            k += 1

            self.data_x.append(f.data_x)
            self.data_y.append(f.data_y)

            ei += f.data_x.size
            self.ranges.append((si,ei))
            si = ei
        self.data_x = numpy.ma.concatenate(self.data_x)
        self.data_y = numpy.ma.concatenate(self.data_y)

    def merge(self, mdlist):
        """
        Merge a list of monodatasets to self
        """
        L = [f for f in self]
        for d in mdlist:
            for f in d:
                L.append(f)
        L = sorted(set(L))
        self.import_framelist(L)

    def pointvalue(self, x):
        """
        Returns a tuple of arrays t,y containing the value of the data
        interpolated at the given point as a function of time.

        We use linear interpolation.
        """
        t = []
        y = []
        for f in self:
            t.append(f.time)
            y.append(f(x))
        return (numpy.array(t), numpy.array(y))

    def purge_nans(self):
        """
        Removes all NaN entries
        """
        L = []
        for f in self:
            f.purge_nans()
            L.append(f)
        self.import_framelist(L)

    def spectrum(self):
        """
        Returns a tuple (x,f,psd) of numpy.arrays

        * x   : 1D array position on the grid
        * f   : 1D array frequency
        * psd : 2D array linear-spectral density of the data at x

        Usage:
          x,f,psd = data.spectrum()
          matplotlib.pyplot.contourf(x,f,psd)

        Note:
          this assumes that the timestep is constant
        """
        xa = self.data_x[self.ranges[0][0]:self.ranges[0][1]].copy()

        dt = self.time[1] - self.time[0]
        df = 1.0 / (dt*(self.nframes-1))
        fa = df * numpy.arange(0, self.nframes)

        psd = []
        for x in xa:
            t, y = self.pointvalue(x)
            (c,s) = numpy.polyfit(t, y, 1)
            y = y - numpy.polyval([c, s], t)
            fft = numpy.fft.fft(y)
            psd.append(numpy.real(fft*fft.conjugate()))

        return (xa, fa, numpy.array(psd))

    def sort(self):
        """
        Sorts the frames within the monodataset
        """
        L = []
        for f in self:
            f.sort()
            L.append(f)
        self.import_framelist(L)

    def time_interp(self, time):
        """
        Linearly interpolate the data to the given time
        """
        idx_a = max(bisect.bisect_right(self.time, time)-1, 0)
        if idx_a == self.nframes - 1:
            return self.frame(self.nframes - 1)
        else:
            idx_b = idx_a + 1

            frame_a = self.frame(idx_a)
            frame_b = self.frame(idx_b)

            tl = 1.0 - (time - self.time[idx_a])/ \
                (self.time[idx_b] - self.time[idx_a])
            frame_i = tl*frame_a + (1.0 - tl)*frame_b

            return frame_i
    def visualize(self):
        """
        Visualize the data using pygraph
        """
        desc,filename = tempfile.mkstemp(suffix=".pyg")
        self.write_pyg(filename)
        os.system("pygraph "+filename)
        os.remove(filename)

    def write(self,filename):
        """
        Writes the data in xgraph format to file
        """
        file = open(filename,"w")
        file.writelines(self.format())
        file.close()
    write_xg = write

    def write_pyg(self, filename):
        """
        Writes the data in pygraph format to file
        """
        dfile = h5py.File(filename, "w")
        dfile['/'].attrs.create('pyg_version', 1, dtype='int32')
        for frm in self:
            dset = dfile.create_dataset(str(frm.index),
                    (2*frm.data_x.shape[0], ))
            x = numpy.array(frm.data_x, dtype='float32')
            y = numpy.array(frm.data_y, dtype='float32')
            dset[0::2] = frm.data_x
            dset[1::2] = frm.data_y
            dset.attrs.create('time', numpy.array((frm.time), dtype='float32'))
        dfile.close()

    def writeframe(self,frame,filename):
        """
        Writes the given frame in gnuplot format
        """
        p = self.frame(frame)
        p.write(filename)

class iterator:
    def __init__(self, data):
        self.frame = -1
        self.data = data

    def __iter__(self):
        return self

    def next(self):
        if self.frame == self.data.nframes-1:
            raise StopIteration
        self.frame += 1
        return self.data.frame(self.frame)
