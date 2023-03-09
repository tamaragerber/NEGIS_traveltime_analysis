# Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream
## Travel-time analysis

These scripts were used to derive the depth-average horizontal anisotropy from airborne radar crosspoints in the Northeast Greenland Ice Stream (NEGIS) which are published in Gerber et al., (2023), Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream, Nature Communications

## Script overview

1) pick reflections in orthogonal wave polarization

	- `pickTWT.m`: script to manually pick travel times for reflections in two polarization directions.
	- crosspoint name and corresponding profiles are read from **crosspoints.xlsx**	
	- save files for each crosspoint with the travel times in **TWTpicks/cp\*_int1_stack0.mat**

2) combine all cp files in TWTpicks to one CP_raw.mat	
	
	- `collect_cp.m` ---> output: **statistics.mat, CP_raw.mat**
	
3) calculate $\Delta \lambda$ with a linear regression.

	- `TWTprocessing.m`: input **CP_raw.mat**

		1) determine reflector depth with DEP-approximated permittivity profile, which is 
		   transformed into velocity

		2) discard reflections where $\Delta$twt is larger than for a monocrystal at the corresponding depth, 
		are shallower than 200 m, or $\Delta$twt is more than 2 standar deviations away from mean trend.

		3) determine depth-average dielectric permittivity for each polarization by linear regression through 0 m, 200 m and 500 m

		4) discard points with less than 5 reflectors.

	- output: **CP.mat, cp_anisotropy.shp**	
	

## Additional data
In order to reproduce the results in the publication above, you additionally need the following data:
- If you want to pick reflections yourself you'll need to download the airborne radar data from PANGAEA (https://doi.org/10.1594/PANGAEA.928569) 
- If you want to continue without picking reflections yourself you can find those on the ERDA data archive (https://doi.org/????) under *raw_data/crosspoint_traveltime_analysis/*
	- the folder **TWT_picks** contains individual files for each crosspoint
	- **CP_raw.mat** contains the raw data of picked reflections
	- **CP.mat** contains the processed data, including the calculated horizontal anisotropy.

Note that you'll need to adjust the lines marked with ***** in the scripts with the corresponding path to the data files. 
