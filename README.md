# NEGIS_traveltime_analysis

These scripts were used to infer horizontal anisotropy from airborne radar crosspoints in the Northeast Greenland Ice Stream (NEGIS) which are published in Gerber et al., (2023), Crystal orientation fabric anisotropy causes directional hardening of the Northeast Greenland Ice Stream, Nature Communications



1) pick reflections in orthogonal polarization

	- pickTWT.m: picks travel times for reflections in two polarization directions.
	- crosspoint name and corresponding profiles are read from crosspoints.xlsx
	- save files for each crosspoint with the travel times in TWTpicks/cp*_int1_stack0.mat

2) combine all cp files in TWTpicks to one CP.mat	
	
	- collect_cp.m ---> output: statistics.mat, CP_raw.mat
	
3) calculate delta lambda with linear regression.

	- TWTprocessing.m: input CP_raw.mat

		1) determine reflector depth with DEP-approximated permittivity profile, which is 
		   transformed into velocity

		2) discard reflections where dtwt is larger than for monocrystal at the corresponding depth 
		   or are shallower than 200 m or dTWT is more than 2*sigma away from mean trend.

		3) determine average epsilon by regression through origin, 200m and 500m

		4) discard points with less than 5 reflectors.

	- output: CP.mat, cp_anisotropy.shp	
	

