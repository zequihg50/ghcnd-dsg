import numpy as np
import netCDF4
import cftime
import math

NC = "ghcn-dsg.nc"
NC = "ghcn-dsg1000.nc"
NC = "ghcn-dsg20240524.nc"

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

class Field:

    def __init__(self, name, nc):
        self.nc = nc
        self.name = name

        self.coordinates = []
        for coord in self.nc.variables[name].getncattr("coordinates").split(" "):
            self.coordinates.append(coord)

        self.axes = {}
        for coord in self.coordinates:
            print(f"Processing axis {coord}.")
            if coord == "time":
                self.axes[coord] = TimeAxisCRA(self.nc)
            else:
                self.axes[coord] = self.nc[coord]

        # the mask is the optimization for subsetting the time values for specific stations
        # see subspace method
        self.mask = np.repeat(True, self.nc.dimensions["timeseries"].size)

    def __str__(self):
        tstart = cftime.num2date(self.axes["time"].min, self.axes["time"].units, self.axes["time"].calendar)
        tend   = cftime.num2date(self.axes["time"].max, self.axes["time"].units, self.axes["time"].calendar)
        return f"{self.name}: Station: {self.axes['station'].shape[0]} Lat: {self.axes['lat'].size}, Lon: {self.axes['lon'].size}, Time: {self.axes['time'].size()} ({tstart}, {tend})\n"

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, item):
        if isinstance(item, slice):
            numdate_start = cftime.date2num(item.start, self.axes["time"].units, self.axes["time"].calendar)
            numdate_stop  = cftime.date2num(item.stop , self.axes["time"].units, self.axes["time"].calendar)

            # create the time slice, since it does not exist in the netCDF
            # assumen time step of 1
            time_slice = np.arange(numdate_start, numdate_stop+1) # includes stop

            # create the numpy array that will hold the output data
            #arr = np.empty((self.nc.dimensions["timeseries"].size, len(time_slice)))
            arr = np.empty((self.mask.sum(), len(time_slice)))

            st = []
            idx = self.axes["time"].idx
            counter = 0 # counts number of stations added, needed to take into account the mask
            for s in idx:
                if not self.mask[s]:
                    continue

                arr[counter] = np.repeat(np.nan, len(time_slice))
                if numdate_start < idx[s][1] or numdate_stop >= idx[s][0]:
                    # locate the values of the field in the netCDF
                    stcsum  = self.nc["rowSize"][...].cumsum()
                    ststart = stcsum[s] - self.nc["rowSize"][s]
                    stend   = stcsum[s]
                    sttimes = self.nc["time"][ststart:stend].astype(np.int32)
                    sttimes_slic = sttimes[(sttimes >= numdate_start) & (sttimes <= numdate_stop)]

                    field_view = self.nc[self.name][ststart:stend]
                    vs = field_view[(sttimes >= numdate_start) & (sttimes <= numdate_stop)]
                    #print(vs.shape) # 6
                    #print(sttimes.shape) # 2500 (el rowSize de la estacion)
                    #print(sttimes_slic.shape) # 6
                    #print(sttimes_slic-numdate_start)
                    arr[counter, sttimes_slic-numdate_start] = vs

                counter += 1

            return arr
        else:
            numdate = cftime.date2num(item, self.axes["time"].units, self.axes["time"].calendar)
            # find the stations with this date
            st = []
            idx = self.axes["time"].idx
            for s in idx:
                if numdate >= idx[s][0] and numdate <= idx[s][1]: # ok, it is in the index
                    # locate the value of the field in the netCDF
                    stcsum  = self.nc["rowSize"][...].cumsum()
                    ststart = stcsum[s] - self.nc["rowSize"][s]
                    stend   = stcsum[s]
                    sttimes = self.nc["time"][ststart:stend]

                    # assumes array is sorted, it should since it is the time coordinate
                    pos = np.searchsorted(sttimes, numdate, side="left")
                    if pos < sttimes.size:
                        if sttimes[pos] == numdate:
                            f = self.nc[self.name]
                            v = f[ststart+pos]
                            st.append((s,v))
                    else:
                        st.append((s, np.nan))

            return np.asarray(st)

    def subspace(self, **kwargs):
        for k in kwargs:
            if k not in self.axes:
                print(f"Invalid axis {k}.")
                return np.array([])

        # if only subspacing time no optimization possible
        if len(kwargs) == 1 and "time" in kwargs:
            return self.__getitem__(kwargs["time"])

        # if subspacing any other coordinates then optimize the time axis
        subset = np.repeat(True, self.nc.dimensions["timeseries"].size)
        for k in kwargs:
            if k == "time":
                continue

            # if subsetting by slice select stations falling in the range
            if k in ["lat", "lon"]:
                if isinstance(kwargs[k], slice):
                    kcoord = self.nc[k]
                    ksubset = (kcoord >= kwargs[k].start) & (kcoord <= kwargs[k].stop)
                else:
                    ksubset = self.nc[k] == kwargs[k]

            if k == "station":
                # ToDo, index by station name
                pass

            subset = np.logical_and(subset, ksubset)

        # now subset contains and array of booleans where True signal stations of interest
        self.mask = subset
        arr = self.__getitem__(kwargs["time"])

        # restore mask and return
        self.mask = np.repeat(True, self.nc.dimensions["timeseries"].size)
        return arr


def read(nc):
    fields = list()
    #for v in ["tasmin", "tasmax", "pr", "snow", "snwd"]:
    for v in ["pr"]:
        fields.append(Field(v, nc))
    return fields

if __name__ == "__main__":
    with netCDF4.Dataset(NC) as nc:
        fields = read(nc)
        print(fields)

        pr = fields[0]
        d1 = cftime.num2date(0, "days since 1916-06-03 00:00:00", "gregorian")
        d2 = cftime.num2date(5, "days since 2018-06-03 00:00:00", "gregorian")
        #print(f"Asking for date: {d1}")
        ##subset = pr.get_time_step(d1)
        #subset = pr[d1]
        #print("Done.")

        #print(f"Asking for slice: {d1}, {d2}")
        #subset = pr[slice(d1,d2)]
        #print(subset)
        #print("Done.")

        #print(f"Subspacing wrong axis.")
        #pr.subspace(oooh=slice(d1,d2))

        lat=np.float32(-29.)
        print(f"Subspacing {d1}-{d2} at lat={lat}.")
        subset = pr.subspace(time=slice(d1,d2), lat=lat) # stations ASN00007125 ASN00007144 in the 1000 dataset (23 stations in the full dataset)
        print(f"Done ({subset.shape}).")

        lat=slice(np.float32(-29.), np.float32(-20.))
        print(f"Subspacing {d1}-{d2} at lat={lat}.")
        subset = pr.subspace(time=slice(d1,d2), lat=lat)
        print(f"Done ({subset.shape}).")

        lat=slice(np.float32(-29.), np.float32(-12.))
        lon=slice(np.float32(117.), np.float32(170.))
        print(f"Subspacing {d1}-{d2} at lat={lat},lon={lon}.")
        subset = pr.subspace(time=slice(d1,d2), lat=lat, lon=lon)
        print(f"Done ({subset.shape}).")
