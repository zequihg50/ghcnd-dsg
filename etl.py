import sys, os
import math, codecs
import numpy as np
import pandas as pd
import netCDF4
import cftime
import requests
from datetime import datetime
from tqdm import tqdm

SRC_DIR = "by_station_3"
CHUNKS = (8192,)

MFLAGS = np.array([ord("A"), ord("B"), ord("D"), ord("H"), ord("K"), ord("L"), ord("O"), ord("P"), ord("T"), ord("W")], dtype="b")
MFLAGS_MEANINGS = " ".join([
"value_in_precipitation_or_snow_is_a_multi-day_total_accumulated_since_last_measurement",
"precipitation_total_formed_from_two_twelve-hour_totals",
"precipitation_total_formed_from_four_six-hour_totals",
"represents_highest_or_lowest_hourly_temperature_or_average_of_hourly_values",
"converted_from_knots",
"temperature_appears_to_be_lagged_with_respect_to_reported_hour_of_observation",
"converted_from_oktas",
"identified_as_missing_presumed_zero_in_DSI_3200_and_3206",
"trace_of_precipitation_snowfall_or_snow_depth",
"converted_from_16-point_WBAN_code_for_wind_direction",
])
QFLAGS = np.array([ord("D"), ord("G"), ord("I"), ord("K"), ord("L"), ord("M"), ord("N"), ord("O"), ord("R"), ord("S"), ord("T"), ord("W"), ord("X"), ord("Z")], dtype="b")
QFLAGS_MEANINGS = " ".join([
"failed_duplicate_check",
"failed_gap_check",
"failed_internal_consistency_check",
"failed_streak_or_frequent-value_check",
"failed_check_on_length_of_multiday_period",
"failed_mega-consistency_check",
"failed_naught_check",
"failed_climatological_outlier_check",
"failed_lagged_range_check",
"failed_spatial_consistency_check",
"failed_temporal_consistency_check",
"temperature_too_warm_for_snow",
"failed_bounds_check",
"flagged_as_a_result_of_an_official_Datzilla_investigation",
])
SFLAGS = np.array([
0, 6, 7, ord("A"), ord("a"), ord("B"), ord("b"), ord("C"), ord("E"), ord("F"), ord("G"), ord("H"), ord("I"), ord("K"),
ord("M"), ord("N"), ord("Q"), ord("R"), ord("r"), ord("S"), ord("s"), ord("T"), ord("U"), ord("u"),ord("W"), ord("X"),
ord("Z"), ord("z")
], dtype="b")
SFLAGS_MEANINGS = " ".join([
"U.S._Cooperative_Summary_of_the_Day_NCDC_DSI-3200",
"CDMP_Cooperative_Summary_of_the_Day_NCDC_DSI-3206",
"U.S._Cooperative_Summary_of_the_Day_Transmitted_via_WxCoder3_NCDC_DSI-3207",
"U.S._Automated_Surface_Observing_System_ASOS_real-time_data_since_2006-01-01",
"Australian_data_from_the_Australian_Bureau_of_Meteorology",
"U.S._ASOS_data_for_October_2000-December_2005_NCDC_DSI-3211",
"Belarus_update",
"Environment_Canada",
"European_Climate_Assessment_and_Dataset_Klein_Tank_et_al_2002",
"U.S._Fort_data",
"Official_Global_Climate_Observing_System_GCOS_or_other_government-supplied_data",
"High_Plains_Regional_Climate_Center_real-time_data",
"International_collection_non_U.S._data_received_through_personal_contacts",
"U.S._Cooperative_Summary_of_the_Day_data_digitized_from_paper_observer_forms_from_2011_to_present",
"Monthly_METAR_Extract_additional_ASOS_data",
"Community_Collaborative_Rain,_Hail,and_Snow_CoCoRaHS",
"Data_from_several_African_countries_that_had_been_quarantined_that_is_withheld_from_public_release_until_permission_was_granted_from_the_respective_meteorological_services",
"NCDC_Reference_Network_Database_Climate_Reference_Network_and_Historical_Climatology_Network-Modernized",
"All-Russian_Research_Institute_of_Hydrometeorological_Information-World_Data_Center",
"Global_Summary_of_the_Day_NCDC_DSI-9618_values_are_derived_from_hourly_synoptic_reports_exchanged_on_the_Global_Telecommunications_System_(GTS)._Daily_values_derived_in_this_fashion_may_differ_significantly_from_true_daily_data,_particularly_for_precipitation._Use_with_caution)",
"China_Meteorological_Administration/National_Meteorological_Information_Center/Climate_Data_Center_cdc.cma.gov.cn",
"SNOwpack_TELemtry_SNOTEL_data_obtained_from_the_Western_Regional_Climate_Center",
"Remote_Automatic_Weather_Station_RAWS_data_obtained_from_the_Western_Regional_Climate_Center",
"Ukraine_update",
"WBAN_ASOS_Summary_of_the_Day_from_NCDC_Integrated_Surface_Data_ISD",
"U.S._First-Order_Summary_of_the_Day_NCDC_DSI-3210",
"Datzilla_official_additions_or_replacements",
"Uzbekistan_update",
])
assert len(MFLAGS) == len(MFLAGS_MEANINGS.split(" "))
assert len(QFLAGS) == len(QFLAGS_MEANINGS.split(" "))
assert len(SFLAGS) == len(SFLAGS_MEANINGS.split(" "))

values_dtype = np.dtype([
    ('station_id', 'S11'),
    ('date', 'S8'),
    ('var', 'S4'),
    ('value', 'i4'),
    ('mflag', 'S1'),
    ('qflag', 'S1'),
    ('sflag', 'S1'),
    ('obstime', 'S4'),
])

MISSING = np.int16(9999)

VS = [
    {
        "name": "PRCP",
        "cfname": "pr",
        "attrs":
            {
                "standard_name": "precipitation_flux",
                "coordinates": "time lat lon alt station",
                "units": "mm",
                "scale_factor": np.float64(.1),
            }
    },
    {
        "name": "TMAX",
        "cfname": "tasmax",
        "attrs":
            {
                "standard_name": "air_temperature",
                "coordinates": "time lat lon alt station",
                "units": "Celsius",
                "scale_factor": np.float64(.1),
            }
    },
    {
        "name": "TMIN",
        "cfname": "tasmin",
        "attrs":
            {
                "standard_name": "air_temperature",
                "coordinates": "time lat lon alt station",
                "units": "Celsius",
                "scale_factor": np.float64(.1),
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
        f.setncattr("featureType", "timeSeries")
        #f.setncattr("cdm_data_type", "timeSeries")
        f.setncattr("Conventions", "CF-1.11")
        
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
            fill_value=np.float64(-999.9),
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
            chunksizes=CHUNKS,  # same chunking than meteorological variables
            compression="zlib",
            complevel=1,
            shuffle=True)
        time.setncattr("standard_name", "time")
        time.setncattr("long_name", "time of measurement")
        time.setncattr("units", "days since 1600-01-01 00:00:00")
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
                chunksizes=CHUNKS,  # by default it chooses 2048, for 1 billion observations involves half a million chunks
                compression="zlib",
                complevel=1,
                shuffle=True,
                fill_value=MISSING)
                #fletcher32=True)  # disabled to avoid issues with kerchunk in the future
            f[v["cfname"]].set_auto_scale(False)
            f[v["cfname"]].setncattr("missing_value", MISSING)
            for attr in v["attrs"]:
                f[v["cfname"]].setncattr(attr, v["attrs"][attr])

            # flags
            mflag = f.createVariable(
                f"{v['cfname']}_mf",
                "b",
                ("obs",),
                chunksizes=CHUNKS,
                compression="zlib",
                complevel=1,
                shuffle=True,
                fill_value=0)
            mflag.setncattr("long_name", "Measurement Flag")
            mflag.setncattr("standard_name", "measurement_flag")
            mflag.setncattr("flag_values", MFLAGS)
            mflag.setncattr("flag_meanings", MFLAGS_MEANINGS)

            qflag = f.createVariable(
                f"{v['cfname']}_qf",
                "b",
                ("obs",),
                chunksizes=CHUNKS,
                compression="zlib",
                complevel=1,
                shuffle=True,
                fill_value=0)
            qflag.setncattr("long_name", "Quality Flag")
            qflag.setncattr("standard_name", "quality_flag")
            qflag.setncattr("flag_values", QFLAGS)
            qflag.setncattr("flag_meanings", QFLAGS_MEANINGS)

            sflag = f.createVariable(
                f"{v['cfname']}_sf",
                "b",
                ("obs",),
                chunksizes=CHUNKS,
                compression="zlib",
                complevel=1,
                shuffle=True,
                fill_value=0)
            sflag.setncattr("long_name", "Source Flag")
            sflag.setncattr("standard_name", "source_flag")
            sflag.setncattr("flag_values", SFLAGS)
            sflag.setncattr("flag_meanings", SFLAGS_MEANINGS)

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
        for i,filename in tqdm(enumerate(os.listdir(SRC_DIR))):
            st = filename.replace(".csv", "")
            f["station"][i] = netCDF4.stringtochar(np.array([st], dtype="S11"))
            lat[i] = stationsdf[stationsdf["station_id"] == st.encode("ascii")].iloc[0]["lat"]
            lon[i] = stationsdf[stationsdf["station_id"] == st.encode("ascii")].iloc[0]["lon"]
            alt[i] = stationsdf[stationsdf["station_id"] == st.encode("ascii")].iloc[0]["elevation"]

            # observations
            data = np.genfromtxt(SRC_DIR + "/" + filename, delimiter=",", dtype=values_dtype)
            all_dates = data["date"]
            uniq_dates = np.sort(np.unique(all_dates)).astype("U8")
            uniq_datetimes = [datetime.strptime(date, "%Y%m%d") for date in uniq_dates]
            uniq_cftimes = [cftime.date2num(date, "days since 1600-01-01 00:00:00", "gregorian") for date in uniq_datetimes]

            frm = int(f["rowSize"][...].sum())
            to = int(frm + len(uniq_cftimes))

            f["time"][frm:to] = uniq_cftimes
            f["rowSize"][i] = len(uniq_cftimes)

            for v in VS:
                vname = v["name"].encode("ascii")
                vdata = data[data["var"] == vname]

                f[v["cfname"]][frm:to] = np.repeat(MISSING, len(uniq_cftimes))
                #f[f"{v['cfname']}_mf"][frm:to] = np.repeat(0, len(uniq_cftimes))

                idx = np.where(
                    np.isin(
                        uniq_dates,
                        vdata["date"].astype("U8"))
                )[0]
                idx = idx + frm
                
                f[v["cfname"]][idx] = vdata["value"]
                f[f"{v['cfname']}_mf"][idx] = vdata["mflag"].view("b")
                f[f"{v['cfname']}_qf"][idx] = vdata["qflag"].view("b")
                f[f"{v['cfname']}_sf"][idx] = vdata["sflag"].view("b")
