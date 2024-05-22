import numpy as np
import netCDF4
import cftime
import math

NC = "ghcn-dsg.nc"
NC = "ghcn-dsg1000.nc"
#NC = "ghcn-dsg20240524.nc"

# Time Axis for Continous Ragged Array
class TimeAxisCRA:

    def __init__(self, nc):
        self.units = nc["time"].getncattr("units")
        self.calendar = nc["time"].getncattr("calendar")

        # the index caches the start and end time points for each station
        self.idx = self.__build_idx(nc)

    def size(self):
        return self.max() - self.min() + 1 # assumes "time step" is 1 (in this case 1 day)

    def min(self):
        v = -1
        for i in self.idx:
            if v > self.idx[i][0] or v == -1:
                v = self.idx[i][0]
        return v

    def max(self):
        v = -1
        for i in self.idx:
            if v < self.idx[i][1] or v == -1:
                v = self.idx[i][1]
        return v

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

        self.axes = {}
        for ax in ["lon", "lat", "alt", "time", "station"]:
            print(f"Processing axis {ax}.")
            if ax == "time":
                self.axes[ax] = TimeAxisCRA(self.nc)
            else:
                self.axes[ax] = self.nc[ax]

class Field:

    def __init__(self, name, nc, dataset):
        self.nc = nc
        self.name = name

        self.axes = {}
        for coord in self.nc.variables[name].getncattr("coordinates").split(" "):
            self.axes[coord] = dataset.axes[coord]

# composition >> inheritance
class FieldView:

    def __init__(self, field, mask=None, time_min=None, time_max=None):
        self.field = field
        self.mask = mask
        self.time_min = time_min
        self.time_max = time_max

        # the mask is the optimization for subsetting the time values for specific stations
        # see subspace method
        if mask is None:
            self.mask = np.repeat(True, self.field.nc.dimensions["timeseries"].size)

        if time_min is None:
            self.time_min = self.field.axes["time"].min()
        if time_max is None:
            self.time_max = self.field.axes["time"].max()

    def __str__(self):
        tstart = cftime.num2date(self.time_min, self.field.axes["time"].units, self.field.axes["time"].calendar)
        tend   = cftime.num2date(self.time_max, self.field.axes["time"].units, self.field.axes["time"].calendar)
        return f"{self.field.name}: Station: {self.mask.sum()} Lat: {self.mask.sum()}, Lon: {self.mask.sum()}, Time: {self.time_max-self.time_min} ({tstart}, {tend})\n"

    def __repr__(self):
        return self.__str__()

    # return the data of the field view as a numpy array
    def data(self):
        # assume time step of 1
        time_slice = np.arange(self.time_min, self.time_max+1) # includes stop

        # create the numpy array that will hold the output data
        arr = np.empty((self.mask.sum(), len(time_slice)))

        st = []
        idx = self.field.axes["time"].idx
        counter = 0 # counts number of stations added, needed to take into account the mask
        for s in idx:
            if not self.mask[s]:
                continue

            arr[counter] = np.repeat(np.nan, len(time_slice))
            if self.time_min < idx[s][1] or self.time_max >= idx[s][0]:
                # locate the values of the field in the netCDF
                stcsum  = self.field.nc["rowSize"][...].cumsum()
                ststart = stcsum[s] - self.field.nc["rowSize"][s]
                stend   = stcsum[s]
                sttimes = self.field.nc["time"][ststart:stend].astype(np.int32)
                sttimes_slic = sttimes[(sttimes >= self.time_min) & (sttimes <= self.time_max)]

                field_view = self.field.nc[self.field.name][ststart:stend]
                vs = field_view[(sttimes >= self.time_min) & (sttimes <= self.time_max)]
                arr[counter, sttimes_slic-self.time_min] = vs

            counter += 1

        return arr

    def subspace(self, **kwargs):
        for k in kwargs:
            if k not in self.field.axes:
                print(f"Invalid axis {k}.")
                return np.array([])

        # if subspacing any other coordinates then optimize the time axis
        #subset = np.repeat(True, self.field.nc.dimensions["timeseries"].size)
        subset = self.mask
        for k in kwargs:
            if k == "time":
                continue

            # if subsetting by slice select stations falling in the range
            if k in ["lat", "lon"]:
                if isinstance(kwargs[k], slice):
                    kcoord = self.field.nc[k]
                    ksubset = (kcoord >= kwargs[k].start) & (kcoord <= kwargs[k].stop)
                else:
                    ksubset = self.field.nc[k] == kwargs[k]

            if k == "station":
                if isinstance(kwargs[k], list):
                    # list of strings, retrieve matches
                    if "_Encoding" in self.field.nc[k].ncattrs():
                        ksubset = np.isin(self.field.nc[k][...], kwargs[k])
                    else:
                        ksubset = np.isin(netCDF4.chartostring(self.field.nc[k][...]), kwargs[k])
                else:
                    # consider a station name only
                    if "_Encoding" in self.field.nc[k].ncattrs():
                        ksubset = self.field.nc[k][...] == kwargs[k]
                    else:
                        ksubset = netCDF4.chartostring(self.field.nc[k][...]) == kwargs[k]

            subset = np.logical_and(subset, ksubset)

        # now subset contains and array of booleans where True signal stations of interest
        new_time_min = self.time_min
        new_time_max = self.time_max
        if "time" in kwargs:
            item = kwargs["time"]
            if not isinstance(item, slice):
                item = slice(item, item)
            numdate_start = cftime.date2num(item.start, self.field.axes["time"].units, self.field.axes["time"].calendar)
            numdate_stop  = cftime.date2num(item.stop , self.field.axes["time"].units, self.field.axes["time"].calendar)
            if numdate_start < self.time_min:
                raise ValueError(f"Asking for data outside the range of the view... View starts at {self.time_min}, asking for {numdate_start}.")
            if numdate_stop > self.time_max:
                raise ValueError(f"Asking for data outside the range of the view... View starts at {self.time_max}, asking for {numdate_stop}.")
            new_time_min = numdate_start
            new_time_max = numdate_stop

        return FieldView(self.field, subset, new_time_min, new_time_max)


def read(nc):
    fields = list()
    dataset = Dataset(nc)
    for v in ["tasmin", "tasmax", "pr", "snow", "snwd"]:
        fields.append(FieldView(Field(v, nc, dataset)))
    return fields

if __name__ == "__main__":
    with netCDF4.Dataset(NC) as nc:
        fields = read(nc)
        print(fields)

        pr = fields[2]
        d1 = cftime.num2date(0, "days since 1916-06-03 00:00:00", "gregorian")
        d2 = cftime.num2date(5, "days since 2018-06-03 00:00:00", "gregorian")

        print(f"Asking for date: {d1}")
        #subset = pr[d1].data()
        subset = pr.subspace(time=d1).data()
        print(f"Done ({subset.shape}).")

        # Too much memory in full dataset. Don't ask for data().
        print(f"Asking for slice: {d1}, {d2}")
        subset = pr.subspace(time=slice(d1,d2))
        print(f"Done ({subset}).")

        print(f"Subspacing wrong axis.")
        pr.subspace(oooh=slice(d1,d2))

        lat=np.float32(-29.)
        print(f"Subspacing {d1}-{d2} at lat={lat}.")
        subset = pr.subspace(time=slice(d1,d2), lat=lat).data() # stations ASN00007125 ASN00007144 in the 1000 dataset (23 stations in the full dataset)
        print(f"Done ({subset.shape}).")

        lat=slice(np.float32(-29.), np.float32(-20.))
        print(f"Subspacing {d1}-{d2} at lat={lat}.")
        subset = pr.subspace(time=slice(d1,d2), lat=lat).data()
        print(f"Done ({subset.shape}).")

        lat=slice(np.float32(-29.), np.float32(-12.))
        lon=slice(np.float32(117.), np.float32(170.))
        print(f"Subspacing {d1}-{d2} at lat={lat},lon={lon}.")
        subset = pr.subspace(time=slice(d1,d2), lat=lat, lon=lon).data()
        print(f"Done ({subset.shape}).")

        stations=['ASN00001028', 'ASN00003002', 'ASN00006049', 'ASN00007064', 'ASN00004064']
        print(f"Subspacing stations {stations}.")
        subset = pr.subspace(station=stations).data()
        print(f"Done ({subset.shape}).")

        print(f"Subspacing stations {stations} at {d1}-{d2}.")
        subset = pr.subspace(station=stations,time=slice(d1,d2)).data()
        print(f"Done ({subset.shape}).")
