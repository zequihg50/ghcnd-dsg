# GHCNd as DSG

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zequihg50/ghcnd-dsg/HEAD?labpath=cf-demo.ipynb)

The Global Historical Climatology Network - Daily (GHCN-Daily) dataset integrates daily climate observations from approximately 30 different data sources.

Download `by_station` data from [here](https://www.ncei.noaa.gov/pub/data/ghcn/daily/). gunzip using `find -type f | parallel gunzip`.

```bash
wget https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt
awk '{print "https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/"$1".csv.gz"}' ghcnd-stations.txt | parallel -j5 wget -nc -q -P by_station/ {}
```

Takes 5 hours in serial for 1 variable and file size is 600MB. Takes 8 hours and 2,5GB for 5 variables.

```
ID = 11 character station identification code
YEAR/MONTH/DAY = 8 character date in YYYYMMDD format (e.g. 19860529 = May 29, 1986)
ELEMENT = 4 character indicator of element type 
DATA VALUE = 5 character data value for ELEMENT 
M-FLAG = 1 character Measurement Flag 
Q-FLAG = 1 character Quality Flag 
S-FLAG = 1 character Source Flag 
OBS-TIME = 4-character time of observation in hour-minute format (i.e. 0700 =7:00 am)
```

## ToDo

- Include flags in the data (1 character type).
- Include metadata from stations (compound data type in station\_info variable?)
- Difference between https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/ and https://www.ncei.noaa.gov/pub/data/ghcn/daily/?
- This is not going to work in the full dataset... See `pr.subspace(T=cf.dt('1916-06-03'))` vs `python field.py`.

In the full dataset:

```
dimensions:
	name_strlen = 11 ;
	timeseries = UNLIMITED ; // (125987 currently)
	obs = UNLIMITED ; // (1122439562 currently)

$ time python field.py 
Processing axis time.
Processing axis lat.
Processing axis lon.
Processing axis alt.
Processing axis station.
[pr: Station: 125987 Lat: 125987, Lon: 125987, Time: 100183 (1750-02-01 00:00:00, 2024-05-17 00:00:00)
]
Subspacing 1916-06-03 00:00:00-2018-06-08 00:00:00 at lat=-29.0.
Done ((23, 37261)).
Subspacing 1916-06-03 00:00:00-2018-06-08 00:00:00 at lat=slice(-29.0, -20.0, None).
Done ((7549, 37261)).
Subspacing 1916-06-03 00:00:00-2018-06-08 00:00:00 at lat=slice(-29.0, -12.0, None),lon=slice(117.0, 170.0, None).
Done ((5367, 37261)).
Subspacing stations ['ASN00001028', 'ASN00003002', 'ASN00006049', 'ASN00007064', 'ASN00004064'].
Done ((5, 100183)).
Subspacing stations ['ASN00001028', 'ASN00003002', 'ASN00006049', 'ASN00007064', 'ASN00004064'] at 1916-06-03 00:00:00-2018-06-08 00:00:00.
Done ((5, 37261)).

real	1m8,879s
user	1m6,967s
sys	0m2,427s
```
