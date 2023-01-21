% script used to pick radar reflections at crosspoints
clear all
close all
path = pwd;

mode = 1;               % 0 = raw data / 1 = spline interpolated 
stacked = 0;            % 0 = single trace / 1 = stacked

Q = readtable('crosspoints.xlsx');
k = 5;                  % choose specific crosspoint

%for cpnum=1:402         % or loop over all crosspoints in file crosspoints.csv

cpnumber = char(Q{k,1});
cpname = char(Q{k,2});
lpname = char(Q{k,3});


% change to radar data directory
cd ****/AWI_radar2018                       % adjust path to radar data 

% load radarfiles
cp = load(cpname);      % flowtransverse
lp = load(lpname);      % flowparallel

% find crosspoint:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lat, lon] = intersections(cp.Latitude,cp.Longitude,lp.Latitude,lp.Longitude);

% troubleshooting: plot radar lines
figure()
plot(cp.Longitude,cp.Latitude);hold on
plot(lp.Longitude,lp.Latitude)
scatter(lon,lat,'ko','filled');
labelpoints(lon,lat,cpnumber,'NE',0.2,1,'color','k','fontsize',10);

if length(lat)>1
    pause()
end

lat = lat(1);
lon = lon(1);
Dcp = distance(cp.Latitude,cp.Longitude,lat,lon,referenceEllipsoid('wgs84'));
[~ ,cpidx] = min(Dcp);

Dlp = distance(lp.Latitude,lp.Longitude,lat,lon,referenceEllipsoid('wgs84'));
[~ ,lpidx] = min(Dlp);

clf

% troubleshooting: radar traces in grayscale
a=subplot('Position',[.1 .1 .4 .8]);
imagesc(lpidx-5:lpidx,lp.Time*1e6,-10*log(lp.Data(:,lpidx-5:lpidx)));hold on
axis([lpidx-5 lpidx lp.Surface(lpidx) lp.Time(end)*1e6]);
title('along flow')
ylabel('TWT [\mu s]')
ya = ylim;

b=subplot('Position',[.5 .1 .4 .8]);
imagesc(cpidx:cpidx+5,cp.Time*1e6,-10*log(cp.Data(:,cpidx:cpidx+5)));hold on
set(b,'Yaxislocation','right')
axis([cpidx cpidx+5 cp.Surface(cpidx) cp.Time(end)*1e6]);
title('across flow')
ylabel('TWT [\mu s]')
yb = ylim;

colormap gray
close all

% convert to power and microseconds
cp_trace = -10*log(cp.Data(:,cpidx));
lp_trace = -10*log(lp.Data(:,lpidx));
cp_time = cp.Time*1e6;
lp_time = lp.Time*1e6;

% set surface and bed
to_cp = cp.Surface(cpidx);
to_lp = lp.Surface(lpidx);
tb_cp = cp.Bottom(cpidx)-to_cp;
tb_lp = lp.Bottom(lpidx)-to_lp;

cd(path)

%%
% interpolation
if mode == 1
   
    cp_time_A = [cp_time(1):0.001:cp_time(end)];
    lp_time_A = [lp_time(1):0.001:lp_time(end)];

    cp_trace_I = interp1(cp_time,double(cp_trace),cp_time_A,'spline');
    lp_trace_I = interp1(lp_time,double(lp_trace),lp_time_A,'spline');

    fig = figure();
    plot(lp_trace,lp_time,'color',[1 0 0 .5]);hold on
    plot(cp_trace,cp_time,'color',[0 0 1 .5]);hold on
    
    cp_trace = cp_trace_I;
    cp_time = cp_time_A;
    lp_trace = lp_trace_I;
    lp_time = lp_time_A;
    
    plot(lp_trace,lp_time,'r:','DisplayName','parallel');hold on
    plot(cp_trace,cp_time,'b:','DisplayName','perpendicular');
    axis ij
    plot([0 500],[tb_lp tb_lp],'r:')
    plot([0 500],[tb_cp tb_cp],'b:')
   
else 
    fig = figure();
    plot(lp_trace,lp_time,'r','DisplayName','parallel');hold on
    plot(cp_trace,cp_time,'b','DisplayName','perpendicular');
end

axis ij
axis([50 450 0 5]);
xlabel('dB')
ylabel('TWT [\mus]')

Help={'Left button: pick line'
    'Right button: delete line'
    'Middle button: switch between cp and lp'
    'Del: Exit'
    's: save'
    'x: adjust axes to other profile'
    'space: reset zoom window'};

% start picking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mode == 0 && stacked ==0
    name = [cpnumber,'_int0_stack0.mat'];
elseif mode == 1 && stacked == 0
    name = [cpnumber,'_int1_stack0.mat'];
elseif mode == 0 && stacked == 1
    name = [cpnumber,'_int0_stack1.mat'];
elseif mode == 1 && stacked == 1
    name = [cpnumber,'_int1_stack1.mat'];
end
%name = [cpnumber,'_twtanalysis.mat'];
try
    load(['TWT_picks/',name]) % load existing reflections if this crosspoint was picked before
    dat = 1;
catch
    dat = 0;
    TWT_cp_bed = [];
    TWT_lp_bed = [];
    x_lp_bed = [];
    x_cp_bed = [];

end

clf
    
if dat ==1
   plot(lp_trace,lp_time-lp_t0,'r','DisplayName','parallel');hold on
   plot(cp_trace,cp_time-cp_t0,'b','DisplayName','perpendicular');
   for i = 1:length(TWT_lp)
    plot(x_lp(i),TWT_lp(i),'ro');hold on
   end
   for i = 1:length(TWT_cp)
    plot(x_cp(i),TWT_cp(i),'bo');hold on
   end
    m = length(TWT_cp);
    n = length(TWT_lp);
else
    lp_t0 = lp.Surface(lpidx)*1e6;
    cp_t0 = cp.Surface(cpidx)*1e6;
    plot(lp_trace,lp_time-lp_t0,'r','DisplayName','parallel');hold on
    plot(cp_trace,cp_time-cp_t0,'b','DisplayName','perpendicular');
end
 

if dat==0
    TWT_cp = [];
    TWT_lp = [];
    n = 0;
    m = 0;
end
    axis ij
    axis([50 450 0 5]);
    xlabel('dB')
    ylabel('TWT [\mus]')

p = 1;

%% start picking reflections

while true
    while p == 1
          set(gcf, 'color', [1 0 0])
          title('parallel to ice flow')
        try 
            [xa,twta,button]=zoomginput(1);
        catch
            bla = 0;
            return
        end
        switch button     
            case 116 % 't' --> pick first arrival
                  
                 [xa,twta,button]=zoomginput(1);
                 var(['h',num2str(n)]) = plot(xa,twta,'r^');hold on
                 pause(1)
                 lp_t0 = lp.Surface(lpidx)*1e6+twta;
                 set(gcf, 'color', [0 0 1])
                 title('perpendicular to ice flow')
                 [xb,twtb,button]=zoomginput(1);
                 var(['h',num2str(n)]) = plot(xb,twtb,'b^');hold on
                 pause(1)
                 cp_t0 =  cp.Surface(cpidx)*1e6+twtb;
                 clf
                 plot(lp_trace,lp_time-lp_t0,'r','DisplayName','parallel');hold on
                 plot(cp_trace,cp_time-cp_t0,'b','DisplayName','perpendicular');
                 plot([0 500],[reflectors_nonan.twt_gprMax_par,reflectors_nonan.twt_gprMax_par],'r')
                 plot([0 500],[reflectors_nonan.twt_gprMax_per,reflectors_nonan.twt_gprMax_per],'b')
                 plot([0 500],[tb_lp tb_lp],'r:')
                 plot([0 500],[tb_cp tb_cp],'b:')
                 axis ij
                 axis([50 450 0 5]);
                 xlabel('dB')
                 ylabel('TWT [\mus]')

            case 32 % 'space' --> reset zoom window
                 y = ylim;
                y1 = y(1);y2 = y(2);
                axis([0 600 y1 y1+1])
               
            case 110 % n
                 n = n+1;
                 [xa,twta,button]=zoomginput(1);
                 var(['h',num2str(n)]) = plot(xa,twta,'xk');hold on
                 pause(1)
                 TWT_lp(n)= NaN;
                 x_lp(n) = NaN;
                 
            case 98 % b --> pick bed reflection
                 [xa,twta,button]=zoomginput(1);
                 var(['h',num2str(n)]) = plot(xa,twta,'r*');hold on
                 pause(1)
                 TWT_lp_bed= twta;
                 x_lp_bed = xa;
                 
            case 117 % u
                y = ylim;
                y1 = y(1);y2 = y(2);
                axis([0 600 y1-(y2-y1) y1]);
                
            case 100 % d
                 y = ylim;
                y1 = y(1);y2 = y(2);
                axis([0 600 y2 y2+(y2-y1)]);
                
            case 1
                n = n+1;
                var(['h',num2str(n)]) = plot(xa,twta,'ro');hold on
                pause(1)
                TWT_lp(n)= twta;
                x_lp(n) = xa;
                
            case 3
                set(var(['h',num2str(n)]),'Visible','off');pause(1)
                TWT_lp(n)= [];
                x_lp(n) = [];
                n = n-1;
                
            case 97
                ya = ylim;
                p = 2;
                break;
                
            case 115 % s
                if n == m
                save(['TWT_picks/',name],...
                    'TWT_cp','x_cp','TWT_lp','x_lp','cp_trace','cp_time','lp_trace',...
                    'lp_time','cp_t0','lp_t0','TWT_cp_bed','TWT_lp_bed','x_cp_bed',...
                    'x_lp_bed','lat','lon','to_cp','to_lp','tb_cp','tb_lp')
                else
                    disp('Warning: dimensions disagree! File is not saved')
                end

            case 127 
                if n == m
                save(['TWT_picks/',name],...
                    'TWT_cp','x_cp','TWT_lp','x_lp','cp_trace','cp_time','lp_trace',...
                    'lp_time','cp_t0','lp_t0','TWT_cp_bed','TWT_lp_bed','x_cp_bed',...
                    'x_lp_bed','lat','lon','to_cp','to_lp','tb_cp','tb_lp')
                else
                    disp('Warning: dimensions disagree! File is not saved')
                end
                p = 3;
                if m==n
                    break;
                end
                
            case 120 % 'x'
                ylim(yb);
                
            case 104
                listdlg('ListString',Help,'SelectionMode','single','Name','Help','ListSize',[300,200]);
        end
    end

    while p == 2
         set(gcf, 'color', [0 0 1])
          title('perpendicular to ice flow')
         try 
            [xb,twtb,button]=zoomginput(1);
         catch
            bla = 0;
            return
         end

        switch button
            case 32 % 'space' --> reset zoom window
                y = ylim;
                y1 = y(1);y2 = y(2);
                axis([0 600 y1 y1+1]);
                
            case 117 % upwards arrow
                  y = ylim;
                y1 = y(1);y2 = y(2);
                axis([0 600 y1-(y2-y1) y1]);
                
            case 110 % n
                 m = m+1;
                 [xb,twtb,button]=zoomginput(1);
                var(['h',num2str(m)]) = plot(xb,twtb,'xk');hold on
                pause(1)
                TWT_cp(m)= NaN;
                x_cp(m) = NaN;
                
            case 98 % b --> pick bed reflection
                 [xb,twtb,button]=zoomginput(1);
                 var(['h',num2str(n)]) = plot(xb,twtb,'b*');hold on
                 pause(1)
                 TWT_cp_bed= twtb;
                 x_cp_bed = xb;
                
            case 100 % downwards arrow
                 y = ylim;
                y1 = y(1);y2 = y(2);
                axis([0 600 y2 y2+(y2-y1)]);
                
            case 1
                m = m+1;
                var(['k',num2str(m)]) =  plot(xb,twtb,'bo');hold on
                pause(1)
                TWT_cp(m)= twtb;
                x_cp(m) = xb;
                
            case 3
                set(var(['k',num2str(m)]),'Visible','off');pause(1)
                TWT_cp(m)= [];
                x_cp(m) = [];
                m = m-1;
                
            case 97
                yb = ylim;
                p = 1;
                break;
                
            case 115
               if n == m
                save(['TWT_picks/',name],...
                    'TWT_cp','x_cp','TWT_lp','x_lp','cp_trace','cp_time','lp_trace',...
                    'lp_time','cp_t0','lp_t0','TWT_cp_bed','TWT_lp_bed','x_cp_bed',...
                    'x_lp_bed','lat','lon','to_cp','to_lp','tb_cp','tb_lp')
                else
                    disp('Warning: dimensions disagree! File is not saved')
                end
            case 127
               if n == m
                save(['TWT_picks/',name],...
                    'TWT_cp','x_cp','TWT_lp','x_lp','cp_trace','cp_time','lp_trace',...
                    'lp_time','cp_t0','lp_t0','TWT_cp_bed','TWT_lp_bed','x_cp_bed',...
                    'x_lp_bed','lat','lon','to_cp','to_lp','tb_cp','tb_lp')
                else
                    disp('Warning: dimensions disagree! File is not saved')
                end
                p = 3;
                if m==n
                    break;
                end
                
            case 120 % 'x'
                ylim(ya);
                
            case 104
                listdlg('ListString',Help,'SelectionMode','single','Name','Help');    
        end
    end
    
    if p == 3
        break
    end
    
end
if n == m
    save(['TWT_picks/',name],...
    'TWT_cp','x_cp','TWT_lp','x_lp','cp_trace','cp_time','lp_trace',...
    'lp_time','cp_t0','lp_t0','TWT_cp_bed','TWT_lp_bed','x_cp_bed',...
    'x_lp_bed','lat','lon','to_cp','to_lp','tb_cp','tb_lp')
else
    disp('Warning: dimensions disagree! File is not saved')
end
% end
