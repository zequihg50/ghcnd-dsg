import numpy as np
import netCDF4
import cftime
import math

NC = "ghcn-dsg.nc"

## Time Axis for Continous Ragged Array
#class TimeAxisCRA:
#
#    def __init__(self, nc):
#        self.units = nc["time"].getncattr("units")
#        self.calendar = nc["time"].getncattr("calendar")
#
#        self.idx = self.__build_idx(nc)
#
#    def __build_idx(self, nc):
#        idx = dict()
#        rs = nc["rowSize"][...]
#        rscsum = nc["rowSize"][...].cumsum()
#
#        chunksize = 1000
#        iters = nc["time"].size / chunksize
#        for i in range(math.ceil(iters)):
#            start, end = i*1000, min((i+1)*1000, nc["time"].size)
#            # for each time point, add the station to the index
#            for j,t in enumerate(nc["time"][start:end]):
#                # find the station of the data point
#                st = 0
#                while start+j >= rscsum[st]:
#                    st+=1
#
#                # add the station to the index
#                if t not in idx:
#                    idx[t] = set()
#                idx[t].add(st)
#
#        return idx

# Time Axis for Continous Ragged Array
class TimeAxisCRA:

    def __init__(self, nc):
        self.units = nc["time"].getncattr("units")
        self.calendar = nc["time"].getncattr("calendar")

        # the index caches the start and end time points for each station
        self.idx = self.__build_idx(nc)

        # cache of the min and max time values
        self.min = -1
        self.max = -1
        for i in self.idx:
            if self.min > self.idx[i][0] or self.min == -1:
                self.min = self.idx[i][0]
            if self.max < self.idx[i][1] or self.min == -1:
                self.max = self.idx[i][1]

    def size(self):
        return self.max - self.min + 1

    def __build_idx(self, nc):
        idx = dict()
        rs = nc["rowSize"][...]
        rscsum = nc["rowSize"][...].cumsum()

        for i in range(nc["station"].shape[0]):
            start = nc["time"][rscsum[i] - rs[i]].item()
            end = nc["time"][rscsum[i] - 1].item()
            idx[i] = (int(start), int(end))

        return idx

class Dataset:

    def __init__(self, nc):
        self.lat = nc["lat"][...]
        self.lon = nc["lon"][...]
        self.stations = nc["station"][...]
        self.time = self.__load_time(nc)

    def __load_time(self, nc):
        return TimeAxisCRA(nc)

    def __str__(self):
        tstart = cftime.num2date(self.time.min, self.time.units, self.time.calendar)
        tend   = cftime.num2date(self.time.max, self.time.units, self.time.calendar)
        return f"Station: {self.stations.shape[0]} Lat: {self.lat.size}, Lon: {self.lon.size}, Time: {self.time.size()} ({tstart}, {tend})\n"

    def __repr__(self):
        return self.__str__()

class Field:

    def __init__(self, dataset, name):
        self.dataset = dataset
        self.name = name

    def __str__(self):
        return f"{self.name}: " + self.dataset.__str__()

    def __repr__(self):
        return self.__str__()

def read(nc):
    dataset = Dataset(nc)
    fields = list()
    for v in ["tasmin", "taxmax"]:
        fields.append(Field(dataset, v))
    return fields

if __name__ == "__main__":
    with netCDF4.Dataset(NC) as nc:
        fields = read(nc)
        print(fields)
