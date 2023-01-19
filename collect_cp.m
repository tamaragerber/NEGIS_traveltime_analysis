% combine picked reflections 

clear all
close all
path = pwd;

egriplat = 75.63;
egriplon = -35.99;
%S = shaperead('../../../../data/AWI_radar2018/egripnor2018_frames_epsg3413.shp');

%%
wgs84 = almanac('earth','wgs84','meters');
cd TWT_picks/
A = dir('cp*_int1*_stack0.mat');
n = length(A);

[cpidx, crossfile, longfile] = textread('../crosspoints.txt','%s %s %s');
cpidx={cpidx'};
k = 0;
cpnames=[];
O={};
k=1;
count = 1;
for i = 1:n
    try
        varName = A(i).name;
        cpname =  ['cp',varName(3:end-16)];
        cpin = find(strcmp([cpidx{:}], cpname)) ;
        cpradar = load(['~/Documents/university_copenhagen/PhD_project/gprMax/gprmax_eastgripmodel/data/AWI_radar2018/',crossfile{cpin}]);
        lpradar = load(['~/Documents/university_copenhagen/PhD_project/gprMax/gprmax_eastgripmodel/data/AWI_radar2018/',longfile{cpin}]);
        [cpxcor cpycor] = polarstereo_fwd(cpradar.Latitude,cpradar.Longitude, 6378137.0,0.08181919,70,-45);
        [lpxcor lpycor] = polarstereo_fwd(lpradar.Latitude,lpradar.Longitude, 6378137.0,0.08181919,70,-45);

        [intx, inty] = intersections(cpxcor,cpycor,lpxcor,lpycor);
        [lat, lon] = polarstereo_inv(intx,inty, 6378137.0,0.08181919,70,-45);
        
        if length(lat)>1
             load(varName,'lat');
             load(varName,'lon');            
             
             if length(lat)>1
                 disp([varName,': check number of crossover']);
                 lat = lat(1);
                 lon= lon(1);
             end
             
             [intx, inty] = polarstereo_fwd(lat,lon, 6378137.0,0.08181919,70,-45);
        end
        Dcp = distance(cpradar.Latitude,cpradar.Longitude,lat(1),lon(1),referenceEllipsoid('wgs84'));
        [~ ,trcidxcp] = min(Dcp);

        Dlp = distance(lpradar.Latitude,lpradar.Longitude,lat(1),lon(1),referenceEllipsoid('wgs84'));
        [~ ,trcidxlp] = min(Dlp);

        dist = distance(cpradar.Latitude(trcidxcp),cpradar.Longitude(trcidxcp),lpradar.Latitude(trcidxlp),lpradar.Longitude(trcidxlp),referenceEllipsoid('wgs84'));
        % find orientation of radar lines
        az_cp = (azimuth(cpradar.Latitude(trcidxcp-1),cpradar.Longitude(trcidxcp-1),cpradar.Latitude(trcidxcp+1),cpradar.Longitude(trcidxcp+1),referenceEllipsoid('wgs84')));
        az_lp = (azimuth(lpradar.Latitude(trcidxlp-1),lpradar.Longitude(trcidxlp-1),lpradar.Latitude(trcidxlp+1),lpradar.Longitude(trcidxlp+1),referenceEllipsoid('wgs84')));
      
         if az_cp>=0 && az_cp<=90
            theta_cp = (90-az_cp);
         elseif az_cp>90 && az_cp<=180
            theta_cp = 270-az_cp;
         elseif az_cp>180 && az_cp<=270
             theta_cp =270-az_cp;
         elseif az_cp>270&& az_cp<=360
             theta_cp= 450-az_cp;
        end
        
         if az_lp>=0 && az_lp<=90
            theta_lp = (90-az_lp);
         elseif az_lp>90 && az_lp<=180
            theta_lp = 270-az_lp;
         elseif az_lp>180 && az_lp<=270
             theta_lp =270-az_lp;
         elseif az_lp>270&& az_lp<=360
             theta_lp= 450-az_lp;
        end

        if any(strcmp(O,['lat',num2str(lat),'lon',num2str(lon)]))
            doubles{1,k} = varName;
            idx = find(strcmp(O,['lat',num2str(lat),'lon',num2str(lon)])==1);
            doubles{2,k} = P{idx};
            k=k+1;
            continue
        else
     
        O{count}= ['lat',num2str(lat),'lon',num2str(lon)];
        P{count} = varName;
        distvec(count) = dist;
        cpnames{count} = cpname;
        CP{count}=load(varName);
        
        CP{count}.theta_lp=theta_lp;
        CP{count}.theta_cp=theta_cp;
    
        CP{count}.angle = abs(theta_cp-theta_lp); % angle between profiles
        
        [CP{count}.TWT_cp,r] = sort(CP{count}.TWT_cp);
        CP{count}.crossflow_datafile = crossfile(cpin);
        CP{count}.alongflow_datafile = longfile(cpin);
        CP{count}.TWT_lp = CP{count}.TWT_lp(r);
        CP{count}.x_cp = CP{count}.x_cp(r);
        CP{count}.x_lp = CP{count}.x_lp(r);

        CP{count}.name = varName;
        CP{count}.cp = ['cp',varName(3:end-16)];
        CP{count}.lat = lat;
        CP{count}.lon = lon;
        
        if length(lat)>1
            stop=0
        end
        
        CP{count}.tracedist = dist;
        CP{count}.trcidxlp = trcidxlp;
        CP{count}.trcidxcp = trcidxcp;
        CP{count}.numberofpicks = length(CP{count}.TWT_lp);
        numberofpicks(count) = length(CP{count}.TWT_lp);
               
        count = count+1;
        end
    catch
        disp([varName,' not found'])
    end    
end
statistics = struct;
statistics.distvec =distvec;
statistics.numberofpicks = numberofpicks;

figure()
subplot(1,2,1)
hist(statistics.distvec,25)
title('distance between traces')
xlabel('[m]')
subplot(1,2,2)
hist(statistics.numberofpicks,25);
title('number of picks')

figure()
scatter(statistics.numberofpicks,statistics.distvec);
xlabel('number of picks')
ylabel('distance between traces')

figure()
for i=1:count-1
    scatter(CP{i}.lon,CP{i}.lat,150,length(CP{i}.TWT_lp),'filled');hold on
    labelpoints(CP{i}.lon,CP{i}.lat,CP{i}.cp,'NE',0.1,1,'color','k','fontsize',8);
end
colorbar()

cd(path)

save('statistics.mat','statistics');
save('CP_raw.mat','CP');

%% save coordinates
for i = 1:count-1
    lat(i) = CP{i}.lat;
    lon(i) =  CP{i}.lon;
    
end
lat = lat';
lon=lon';

coord = [lat,lon];
save('cp_coordinates.mat','coord')