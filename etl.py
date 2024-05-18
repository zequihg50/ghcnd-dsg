import sys, os
import math, codecs
import numpy as np
import pandas as pd
import netCDF4
import cftime
import requests
from datetime import datetime
from tqdm import tqdm

values_dtype = np.dtype([
    ('station_id', 'S11'),
    ('date', 'S8'),
    ('var', 'S4'),
    ('value', 'f4'),
    ('mflag', 'S1'),
    ('qflag', 'S1'),
    ('sflag', 'S1'),
    ('obstime', 'S4'),
])

MISSING = np.float32(999.9)

VS = [
    {
        "name": "PRCP",
        "cfname": "pr",
        "attrs":
            {
                "standard_name": "precipitation_flux",
                "coordinates": "time lat lon alt station",
                "units": "kg m-2 s-1"
            }
    },
    {
        "name": "TMAX",
        "cfname": "tasmax",
        "attrs":
            {
                "standard_name": "air_temperature",
                "coordinates": "time lat lon alt station",
                "units": "K",
            }
    },
    {
        "name": "TMIN",
        "cfname": "tasmin",
        "attrs":
            {
                "standard_name": "air_temperature",
                "coordinates": "time lat lon alt station",
                "units": "K",
            }
    },
    {
        "name": "SNOW",
        "cfname": "snow",
        "attrs":
            {
                "standard_name": "snow_amount",
                "coordinates": "time lat lon alt station",
                "units": "mm",
            }
    },
    {
        "name": "SNWD",
        "cfname": "snwd",
        "attrs":
            {
                "standard_name": "snow_depth",
                "coordinates": "time lat lon alt station",
                "units": "mm",
            }
    },
]

def parse_station(line):
    station = tuple([
        line[0:11].rstrip(),
        float(line[12:20].rstrip()),
        float(line[21:30].rstrip()),
        float(line[31:37].rstrip()),
        line[38:40].rstrip(),
        line[41:71].rstrip(),
        line[72:75].rstrip(),
        line[76:79].rstrip(),
        line[80:85].rstrip('\n').rstrip(),
    ])

    return station

def parse_stations(f):
    for line in codecs.open(f, 'r', encoding='utf-8', errors='ignore'):
        station = parse_station(line)
        yield station

def parse_date(d):
    return datetime.date(
        int(d[0:4]),
        int(d[4:6]),
        int(d[6:8]))
        
if __name__ == "__main__":
    with netCDF4.Dataset("ghcn-dsg.nc", "w") as f:
        f.setncattr("featureType", "TimeSeries")
        f.setncattr("cdm_data_type", "TimeSeries")
        f.setncattr("Conventions", "COARDS, CF-1.6, ACDD-1.3")
        
        f.createDimension("name_strlen", 11)
        f.createDimension("timeseries", None)

        lon = f.createVariable(
            "lon",
            "f4",
            ("timeseries",),
            compression="zlib",
            complevel=1,
            shuffle=True)
        lon.setncattr("standard_name", "longitude")
        lon.setncattr("long_name", "station longitude")
        lon.setncattr("units", "degrees_east")
        lon.setncattr("axis", "X")
        lon.setncattr("_CoordinateAxisType", "Longitude")
        
        lat = f.createVariable(
            "lat",
            "f4",
            ("timeseries",),
            compression="zlib",
            complevel=1,
            shuffle=True)
        lat.setncattr("standard_name", "latitude")
        lat.setncattr("long_name", "station latitude")
        lat.setncattr("units", "degrees_north")
        lat.setncattr("axis", "Y")
        lat.setncattr("_CoordinateAxisType", "Latitude")
            
        alt = f.createVariable(
            "alt",
            "f4",
            ("timeseries",),
            compression="zlib",
            complevel=1,
            shuffle=True)
        alt.setncattr("long_name", "vertical distance above the surface")
        alt.setncattr("standard_name", "height")
        alt.setncattr("units", "m")
        alt.setncattr("positive", "up")
        alt.setncattr("axis", "Z")
        alt.setncattr("_CoordinateAxisType", "Height")
        
        station = f.createVariable(
            "station",
            "S1",
            ("timeseries","name_strlen"),
            compression="zlib",
            complevel=1,
            shuffle=True)
        station.setncattr("long_name", "station name")
        station.setncattr("long_name", "Station")
        station.setncattr("cf_role", "timeseries_id")

        obs = f.createDimension("obs", None)

        time = f.createVariable(
            "time",
            "f4",
            ("obs",),
            compression="zlib",
            complevel=1,
            shuffle=True)
        time.setncattr("standard_name", "time")
        time.setncattr("long_name", "time of measurement")
        time.setncattr("units", "days since 1500-01-01 00:00:00")
        time.setncattr("calendar", "gregorian")
        time.setncattr("axis", "T")
        time.setncattr("_CoordinateAxisType", "Time")

        rowSize = f.createVariable(
            "rowSize",
            "i4",
            ("timeseries",),
            compression="zlib",
            complevel=1,
            shuffle=True)
        rowSize.setncattr("long_name", "Number of Observations for this TimeSeries")
        rowSize.setncattr("sample_dimension", "obs")
        
        for v in VS:
            f.createVariable(
                v["cfname"],
                MISSING.dtype,
                ("obs",),
                compression="zlib",
                complevel=1,
                shuffle=True,
                fill_value=MISSING,
                fletcher32=True)
            f[v["cfname"]].setncattr("missing_value", MISSING)
            for attr in v["attrs"]:
                f[v["cfname"]].setncattr(attr, v["attrs"][attr])

        # stations
        stations_dtype = np.dtype([
                    ('station_id', 'S11'),
                    ('lat', 'f8'),
                    ('lon', 'f8'),
                    ('elevation', 'f4'),
                    ('state', 'S2'),
                    #('name', h5py.string_dtype("utf-8", 30)),
                    ('name', 'U30'),
                    ('gsn_flag', 'S3'),
                    ('hcn_crn_flag', 'S3'),
                    ('wmo_id', 'S5'),
        ])
        stations = np.array([x for x in parse_stations("ghcnd-stations.txt")], dtype=stations_dtype)
        stationsdf = pd.DataFrame(stations)

        # write
        f["rowSize"][:] = 0
        for i,filename in tqdm(enumerate(os.listdir("by_station"))):
            st = filename.replace(".csv", "")
            f["station"][i] = netCDF4.stringtochar(np.array([st], dtype="S11"))
            lat[i] = stationsdf[stationsdf["station_id"] == st.encode("ascii")].iloc[0]["lat"]
            lon[i] = stationsdf[stationsdf["station_id"] == st.encode("ascii")].iloc[0]["lon"]
            alt[i] = stationsdf[stationsdf["station_id"] == st.encode("ascii")].iloc[0]["elevation"]

            # observations
            data = np.genfromtxt("by_station/" + filename, delimiter=",", dtype=values_dtype)
            all_dates = data["date"]
            uniq_dates = np.sort(np.unique(all_dates)).astype("U8")
            uniq_datetimes = [datetime.strptime(date, "%Y%m%d") for date in uniq_dates]
            uniq_cftimes = [cftime.date2num(date, "days since 1500-01-01 00:00:00", "gregorian") for date in uniq_datetimes]

            frm = int(f["rowSize"][...].sum())
            to = int(frm + len(uniq_cftimes))

            f["time"][frm:to] = uniq_cftimes
            f["rowSize"][i] = len(uniq_cftimes)

            for v in VS:
                vname = v["name"].encode("ascii")
                vdata = data[data["var"] == vname]

                f[v["cfname"]][frm:to] = np.repeat(MISSING, len(uniq_cftimes))
                
                idx = np.where(
                    np.isin(
                        uniq_dates,
                        vdata["date"].astype("U8"))
                )[0]
                idx = idx + frm
                
                f[v["cfname"]][idx] = vdata["value"]
