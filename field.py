import numpy as np
import netCDF4
import math

NC = "ghcn-dsg.nc"

# Time Axis for Continous Ragged Array
class TimeAxisCRA:

    def __init__(self, nc):
        self.units = nc["time"].getncattr("units")
        self.calendar = nc["time"].getncattr("calendar")

        self.idx = self.__build_idx(nc)

    def __build_idx(self, nc):
        idx = dict()
        rs = nc["rowSize"][...]
        rscsum = nc["rowSize"][...].cumsum()

        chunksize = 1000
        iters = nc["time"].size / chunksize
        for i in range(math.ceil(iters)):
            start, end = i*1000, min((i+1)*1000, nc["time"].size)
            # for each time point, add the station to the index
            for j,t in enumerate(nc["time"][start:end]):
                # find the station of the data point
                st = 0
                while start+j >= rscsum[st]:
                    st+=1

                # add the station to the index
                if t not in idx:
                    idx[t] = set()
                idx[t].add(st)

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
        return f"Station: {self.stations.shape[0]} Lat: {self.lat.size}, Lon: {self.lon.size}, Time: {len(self.time.idx)}\n"

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
