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
