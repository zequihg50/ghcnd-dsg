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
        return self.max - self.min + 1 # assumes "time step" is 1 (in this case 1 day)

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
        self.nc = nc
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

    def nc(self):
        return self.dataset.nc

    def get_time_step(self, date):
        numdate = cftime.date2num(date, self.dataset.time.units, self.dataset.time.calendar)
        # find the stations with this date
        st = []
        idx = self.dataset.time.idx
        for s in idx:
            if numdate >= idx[s][0] and numdate <= idx[s][1]: # ok, it is in the index
                # locate the value of the field in the netCDF
                stcsum  = self.nc()["rowSize"][...].cumsum()
                ststart = stcsum[s] - self.nc()["rowSize"][s]
                stend   = stcsum[s]
                sttimes = self.nc()["time"][ststart:stend]

                # assumes array is sorted, it should since it is the time coordinate
                pos = np.searchsorted(sttimes, numdate, side="left")
                if pos < sttimes.size:
                    if sttimes[pos] == numdate:
                        f = self.nc()[self.name]
                        v = f[ststart+pos]
                        st.append((s,v))
                else:
                    st.append((s, np.nan))

        return np.asarray(st)

    def __getitem__(self, item):
        if isinstance(item, slice):
            numdate_start = cftime.date2num(item.start, self.dataset.time.units, self.dataset.time.calendar)
            numdate_stop  = cftime.date2num(item.stop , self.dataset.time.units, self.dataset.time.calendar)

            # create the time slice, since it does not exist in the netCDF
            # assumen time step of 1
            time_slice = np.arange(numdate_start, numdate_stop+1) # includes stop

            # create the numpy array that will hold the output data
            arr = np.empty((self.nc().dimensions["timeseries"].size, len(time_slice)))

            st = []
            idx = self.dataset.time.idx
            for s in idx:
                arr[s] = np.repeat(np.nan, len(time_slice))
                if numdate_start < idx[s][1] or numdate_stop >= idx[s][0]:
                    # locate the values of the field in the netCDF
                    stcsum  = self.nc()["rowSize"][...].cumsum()
                    ststart = stcsum[s] - self.nc()["rowSize"][s]
                    stend   = stcsum[s]
                    sttimes = self.nc()["time"][ststart:stend].astype(np.int32)
                    sttimes_slic = sttimes[(sttimes >= numdate_start) & (sttimes <= numdate_stop)]

                    field_view = self.nc()[self.name][ststart:stend]
                    vs = field_view[(sttimes >= numdate_start) & (sttimes <= numdate_stop)]
                    #print(vs.shape) # 6
                    #print(sttimes.shape) # 2500 (el rowSize de la estacion)
                    #print(sttimes_slic.shape) # 6
                    #print(sttimes_slic-numdate_start)
                    arr[s, sttimes_slic-numdate_start] = vs

            return arr

    def __str__(self):
        return f"{self.name}: " + self.dataset.__str__()

    def __repr__(self):
        return self.__str__()

def read(nc):
    dataset = Dataset(nc)
    fields = list()
    for v in ["tasmin", "taxmax", "pr"]:
        fields.append(Field(dataset, v))
    return fields

if __name__ == "__main__":
    with netCDF4.Dataset(NC) as nc:
        fields = read(nc)
        print(fields)

        pr = fields[2]
        d1 = cftime.num2date(0, "days since 1916-06-03 00:00:00", "gregorian")
        d2 = cftime.num2date(5, "days since 1916-06-03 00:00:00", "gregorian")
        print(f"Asking for date: {d1}")
        subset = pr.get_time_step(d1)
        print("Done.")

        print(f"Asking for slice: {d1}, {d2}")
        subset = pr[slice(d1,d2)]
        print(subset)
        print("Done.")
