clear all
close all
warning off
% Initializing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

egriplat = 75.63;
egriplon = -35.99;

load('CP_raw.mat');

%%
% TWT analysis: reflector depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n= length(CP);

% determine dTWT
% .........................................................................
for i = 1:n
    CP{i}.dTWT = CP{i}.TWT_cp-CP{i}.TWT_lp;
end


% travel time - pre-processing:
% .........................................................................
load 'permittivity'                                                     % permittivity profile at EGRIP
v = (physconst('lightspeed')./sqrt(permittivity.perm)).*1e-6;           % assumed velocity profile 

TWTm(1) = 2*permittivity.depth(1)/v(1);                                 
 for i = 2:length(v)
     TWTm(i) = TWTm(i-1)+2*(permittivity.depth(i)-permittivity.depth(i-1))/v(i);
 end
 
TWTm = TWTm';
depth = [0:0.1:3500];                                                           % depth
TWT = interp1([0;permittivity.depth],[0;TWTm],depth,'linear','extrap');         % travel time, in micros

% determine max dTWT
% .........................................................................
TWTmax = 2*depth*0.0096*1e6/physconst('lightspeed');                             % maximum dTWT in micros. 

% discard double points
% .........................................................................
for i=1:n
   L(i,1) = CP{i}.lat;
   L(i,2) = CP{i}.lon;  
end
[B, iB, iA] = unique(L,'rows');


for i = 1:n
    try
       CP{i}.meanTWT = mean([CP{i}.TWT_lp;CP{i}.TWT_cp]);                           % average between parallel and perpendicular trace     
       for j = 1:length(CP{i}.meanTWT)
            [~,idx] = min(abs(CP{i}.meanTWT(j)-TWT));                               % find closest depth
            [~,idx2] = min(abs(CP{i}.TWT_lp(j)-TWT));  
            [~,idx3] = min(abs(CP{i}.TWT_cp(j)-TWT));  
            twt_diff = abs(CP{i}.meanTWT(j)-TWT(idx));                              % check for errors in twt-depth conversion
            twt_diff2 = abs(CP{i}.TWT_lp(j)-TWT(idx2)); 
            twt_diff3 = abs(CP{i}.TWT_cp(j)-TWT(idx3)); 
            
            % discard picked reflectors where dTWT is larger than
            % monocrystal
            %..............................................................
            if abs(CP{i}.dTWT(j))>TWTmax(idx)                                       
                d(j) = NaN;
                d_lp(j) = NaN;
                d_cp(j) = NaN;
            else
                % estimate depth
                % ........................................................
                d(j) = depth(idx);
                d_lp(j) = depth(idx2);
                d_cp(j) = depth(idx3);
            end
        end
    CP{i}.depth = d;
    CP{i}.depth_lp = d_lp;
    CP{i}.depth_cp = d_cp;
    clear d kk qq d_lp d_cp
  
    % discard picks where dTWT > dTWTmax, as they are probably wrong
    %......................................................................
    kk = find(isnan(CP{i}.depth)==1);
    CP{i}.meanTWT(kk) = [];
    CP{i}.TWT_cp(kk) = [];
    CP{i}.TWT_lp(kk) = [];
    CP{i}.x_cp(kk) = [];
    CP{i}.x_lp(kk) = [];
    CP{i}.dTWT(kk) = [];
    CP{i}.depth(kk) = [];
    CP{i}.depth_cp(kk) = [];
    CP{i}.depth_lp(kk) = [];
    
    % discard picks shallower than 200m, since sampling interval is
    % larger than expected dTWT
    % .....................................................................
    qq = find((CP{i}.depth)<200);
    CP{i}.meanTWT(qq) = [];
    CP{i}.TWT_cp(qq) = [];
    CP{i}.TWT_lp(qq) = [];
    CP{i}.x_cp(qq) = [];
    CP{i}.x_lp(qq) = [];
    CP{i}.dTWT(qq) = [];
    CP{i}.depth(qq) = [];
    CP{i}.depth_cp(qq) = [];
    CP{i}.depth_lp(qq) = [];
        
    % remove outliers: picks where dTWT is more than 2sigma away from mean
    % trend
    % .....................................................................
    det = detrend(CP{i}.dTWT);
    s = std(det,'omitnan');
    m = mean(det,'omitnan');
    outlier = find(abs(det)>m+2*s);
    CP{i}.meanTWT(outlier) = [];
    CP{i}.TWT_cp(outlier) = [];
    CP{i}.TWT_lp(outlier) = [];
    CP{i}.x_cp(outlier) = [];
    CP{i}.x_lp(outlier) = [];
    CP{i}.dTWT(outlier) = [];
    CP{i}.depth(outlier) = [];
    CP{i}.depth_cp(outlier) = [];
    CP{i}.depth_lp(outlier) = [];
    det = detrend(CP{i}.dTWT);
    
    clear outlier
    
    catch
        disp(['catched ', CP{i}.cp])
    end
end

%% check depth differences:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure();
% for i= 1:n
%     CP{i}.depth_diff = CP{i}.depth_cp-CP{i}.depth_lp;
%     M = max(CP{i}.depth_diff);
%     if isempty(M)==1
%         M=NaN;
%     end
%     scatter(CP{i}.lon,CP{i}.lat,25,M,'filled');hold on
% end
% colorbar()
% title('maximum depth difference between cp and lp [m]')


%% find 6 cp close to EastGRIP
%..........................................................................
wgs84 = almanac('earth','wgs84','meters');

for i = 1:n
    egrip_distance(i) = distance(CP{i}.lat,CP{i}.lon,egriplat,egriplon,wgs84);
    CP{i}.egrip_distance = egrip_distance(i);
end
[a, b]= sort(egrip_distance);
egriptraces_idx = b(1:6); 

%% determine average epsilon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i =1:n
    % linear regression through origin
    %......................................................................
    B = [0;CP{i}.depth(:)]\[0;CP{i}.TWT_cp(:)];                             % regression through origin; B= slope
    CP{i}.TWT_cps = [CP{i}.depth(:)]*B;                                     % Calculate Fitted Line
    CP{i}.slope_eps_cp = (B*physconst('lightspeed')*1e-6/2)^2;              % depth-average epsilon across-flow    
    
    B = [0;CP{i}.depth(:)]\[0;CP{i}.TWT_lp(:)];                             % regression through origin
    CP{i}.TWT_lps = [CP{i}.depth(:)]*B;                                     % Calculate Fitted Line
    CP{i}.slope_eps_lp = (B*physconst('lightspeed')*1e-6/2)^2;              % depth-average epsilon along-flow
    
    CP{i}.slope_deps = abs(CP{i}.slope_eps_cp-CP{i}.slope_eps_lp);          % depth-average difference in permittivity
    CP{i}.slope_dlambda = CP{i}.slope_deps/0.034;                           % depth-average horizontal eigenvalue difference
    CP{i}.dTWTs = CP{i}.TWT_cps-CP{i}.TWT_lps;
    
    % linear regression through 200m point
    % .....................................................................
    CP{i}.t0_200_cp =  sqrt(CP{i}.slope_eps_cp)*400/(physconst('lightspeed')*1e-6);
    CP{i}.t0_200_lp =  sqrt(CP{i}.slope_eps_lp)*400/(physconst('lightspeed')*1e-6);
    
    q = find(CP{i}.depth>200);
    B = (CP{i}.depth(q)-200)'\(CP{i}.TWT_cp(q)-CP{i}.t0_200_cp)';           % regression through 200m
    CP{i}.TWT_cps200 = [CP{i}.depth(q)]*B;                                  % Calculate Fitted Line
    CP{i}.slope_eps_cp200 = (B*physconst('lightspeed')*1e-6/2)^2;
    
    B = (CP{i}.depth(q)-200)'\(CP{i}.TWT_lp(q)-CP{i}.t0_200_lp)';           % regression through 200m
    CP{i}.TWT_lps200 = [CP{i}.depth(q)]*B;                                  % Calculate Fitted Line
    CP{i}.slope_eps_lp200 = (B*physconst('lightspeed')*1e-6/2)^2;
    CP{i}.depth_200m = CP{i}.depth(q); 

    CP{i}.slope_deps200 = abs(CP{i}.slope_eps_cp200-CP{i}.slope_eps_lp200);
    CP{i}.slope_dlambda200 = CP{i}.slope_deps200/0.034;
    
    % linear regression through 500m point
    % .....................................................................
    CP{i}.t0_500_cp =  sqrt(CP{i}.slope_eps_cp)*1000/(physconst('lightspeed')*1e-6);
    CP{i}.t0_500_lp =  sqrt(CP{i}.slope_eps_lp)*1000/(physconst('lightspeed')*1e-6);   
    
    q = find(CP{i}.depth>500);
    B = (CP{i}.depth(q)-500)'\(CP{i}.TWT_cp(q)-CP{i}.t0_500_cp)';           % regression through origin
    CP{i}.TWT_cps500 = [CP{i}.depth(q)]*B;                                  % Calculate Fitted Line
    CP{i}.slope_eps_cp500 = (B*physconst('lightspeed')*1e-6/2)^2;
    
    B = (CP{i}.depth(q)-500)'\(CP{i}.TWT_lp(q)-CP{i}.t0_500_lp)';           % regression through origin
    CP{i}.TWT_lps500 = [CP{i}.depth(q)]*B;                                  % Calculate Fitted Line
    CP{i}.slope_eps_lp500 = (B*physconst('lightspeed')*1e-6/2)^2;
    CP{i}.depth_500m = CP{i}.depth(q); 
    
    CP{i}.slope_deps500 = abs(CP{i}.slope_eps_cp500-CP{i}.slope_eps_lp500);
    CP{i}.slope_dlambda500 = CP{i}.slope_deps500/0.034;
    
    % linear regression NOT through origin
    % .....................................................................
    p1 = polyfit([0,CP{i}.depth],[0,CP{i}.TWT_cp],1);
    CP{i}.TWT_cps_y = polyval(p1,[CP{i}.depth]);                            % regression not going through origin
    CP{i}.slope_eps_cp_y = (p1(1)*physconst('lightspeed')*1e-6/2)^2;
     
    p2 =polyfit([0,CP{i}.depth],[0,CP{i}.TWT_lp],1);
    CP{i}.TWT_lps_y = polyval(p2,[CP{i}.depth]);
    CP{i}.slope_eps_lp_y = (p2(1)*physconst('lightspeed')*1e-6/2)^2;
    
    CP{i}.slope_deps_y = abs(CP{i}.slope_eps_cp_y-CP{i}.slope_eps_lp_y);
    CP{i}.slope_dlambda_y = CP{i}.slope_deps_y/0.034;
    
    % calculate depth-varying epsilon
    % .....................................................................
    CP{i}.eps_cpav = (physconst('lightspeed')*1e-6*CP{i}.TWT_cp./(2*CP{i}.depth)).^2;
    CP{i}.eps_lpav = (physconst('lightspeed')*1e-6*CP{i}.TWT_lp./(2*CP{i}.depth)).^2;
    CP{i}.deps_av = CP{i}.eps_cpav-CP{i}.eps_lpav;

    CP{i}.dlambda_av = CP{i}.deps_av/0.034;
    CP{i}.depth_av = CP{i}.depth*0.5;
              
    % raw picks
    %......................................................................
    for k = 2:length(CP{i}.TWT_cp)
        % raw picks
        %..........................................................
        CP{i}.eps_cp(k-1) = (physconst('lightspeed')*1e-6*(CP{i}.TWT_cp(k)-CP{i}.TWT_cp(k-1))./(2*(CP{i}.depth(k)-CP{i}.depth(k-1)))).^2;
        CP{i}.eps_lp(k-1) = (physconst('lightspeed')*1e-6*(CP{i}.TWT_lp(k)-CP{i}.TWT_lp(k-1))./(2*(CP{i}.depth(k)-CP{i}.depth(k-1)))).^2;
    end
    try
        CP{i}.deps = CP{i}.eps_cp - CP{i}.eps_lp;
        CP{i}.dlambda = CP{i}.deps/0.034;
        CP{i}.dlambda_rawav = abs(mean(CP{i}.dlambda,'omitnan'));
    catch
        CP{i}.deps = NaN;
        CP{i}.dlambda = NaN;
        CP{i}.dlambda_rawav = NaN;
    end
    
    % discard points with less than 5 picked reflectors
    %......................................................................
    if CP{i}.numberofpicks<5
        CP{i}.deps=NaN;
        CP{i}.dlambda = NaN;
        CP{i}.slope_dlambda = NaN;
        CP{i}.dlambda_rawav =NaN;
    end
end

%% PLOTTING TWT-difference

c=colormap(lines);

figure()
xlabel('\Delta t')
ylabel('time [\mu s]')

j=0;
l=1;
 for i = egriptraces_idx
     j=j+1;
     plot(CP{i}.dTWT,CP{i}.TWT_cp,':*','color',c(j,:));hold on
     disp(CP{i}.name)
     l=l+1;
 end

legend('cp217','cp421','cp419','cp138','cp136','cp420')
set(gca,'fontsize',14)
xlabel('\Delta TWT [\mus]')
ylabel('TWT_cp [\mu s]')
title('\Delta TWT')
axis ij 



%% Plotting: 0-regression vs. regression trough 200/500m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
figure()
subplot(1,2,1)
j=0;
colormap lines
plot(egriplon,egriplat,'k*');hold on
for i = egriptraces_idx
    j=j+1;
    scatter(CP{i}.lon,CP{i}.lat,25,c(j,:),'filled');hold on
end
axis square
box on
xlabel('Longitude')
ylabel('Latitude')
legend('EastGRIP');
title('Crosspoint location')
grid on

subplot(1,2,2)

j=0;
for i = egriptraces_idx
     j=j+1;
     plot(CP{i}.dTWT,CP{i}.depth,':*','color',c(j,:));hold on
     plot([0;CP{i}.dTWTs],[0,CP{i}.depth],'-','color',c(j,:),'linewidth',2);hold on
     plot(CP{i}.TWT_cps200 -CP{i}.TWT_lps200,[CP{i}.depth],':','color',c(j,:),'linewidth',2);  
     plot(CP{i}.TWT_cps500 -CP{i}.TWT_lps500,[CP{i}.depth_500m],'--','color',c(j,:),'linewidth',2);  
 end

legend('CP \Delta TWT picks','CP 0-regression','CP 200m regression','CP 500m regression')
xlabel('\Delta TWT [\mus]')
ylabel('depth [m]')
title('\Delta TWT')
axis ij
axis([-0.02 0.1 0 3000])

%% Create table with CP close to EGRIP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 j=0;
 for i = egriptraces_idx
     j=j+1;
     disp(CP{i}.name)
     tab(j,1) = CP{i}.lat;
     tab(j,2) = CP{i}.lon;
     tab(j,3) = CP{i}.slope_dlambda;
     tab(j,4) = CP{i}.slope_dlambda200;
     tab(j,5) = CP{i}.slope_dlambda500;
     if CP{i}.angle>90
         tab(j,6) =180- CP{i}.angle;
     else
        tab(j,6) = CP{i}.angle;
     end
 end

%% save output
%..........................................................................
save('CP.mat','CP')

%% create shapefile
%..........................................................................
Data.Geometry = 'Polygon' ;
m=1;
wgs84 = almanac('earth','wgs84','meters');
clear L lat lon X Y V W M label J K D
for i = 1:n
    try
       if length(CP{i}.TWT_cp) > 5
            ble=0; 
            [x y] = polarstereo_fwd(CP{i}.lat,CP{i}.lon, 6378137.0,0.08181919,70,-45);
            X(i) = x;    
            Y(i) = y ;  
            lat(i) = CP{i}.lat;
            lon(i) = CP{i}.lon;
            V(i) = abs(CP{i}.slope_dlambda(1));
            K(i) = CP{i}.slope_dlambda200(1);
            A(i) = CP{i}.slope_dlambda500(1);
            J(i) = CP{i}.dlambda_rawav;
            O(i) = CP{i}.slope_dlambda_y(1);
            D(i) = CP{i}.depth(end);
            
            W(i) = CP{i}.angle;
            L(i) = CP{i}.theta_lp;
            M(i) = CP{i}.theta_cp;
            label(i)= str2double(CP{i}.cp(3:end));
            num(i) = CP{i}.numberofpicks;
            m= m+1;
       else
           [x y] = polarstereo_fwd(CP{i}.lat,CP{i}.lon, 6378137.0,0.08181919,70,-45);
           X(i) = x;    % latitude
           Y(i) = y ;  % longitude
           lat(i) = CP{i}.lat;
           lon(i) = CP{i}.lon;
           V(i) = -100;
           K(i) = -100;
           A(i) = -100;
           J(i) = -100;
           O(i) = -100;
            W(i) = CP{i}.angle;
            D(i) = CP{i}.depth(end);
          L(i) = CP{i}.theta_lp;
          M(i) = CP{i}.theta_cp;
          label(i)= str2double(CP{i}.cp(3:end));
          num(i) = CP{i}.numberofpicks;
       end
    catch
        disp([CP{i}.name,' not found']);
         
    end
end

p = mappoint(X,Y,'dl_0',V,'dl_200',K,'dl_500',A,'dl_rawav',J,'dl_unc',O,'angle',W,'theta_lp',L,'theta_cp',M,'name',label,'picks',num,'depth',D);
Data.Name = 'cp' ;   % some attribute/ name
shapewrite(p, 'cp_anisotropy.shp')
