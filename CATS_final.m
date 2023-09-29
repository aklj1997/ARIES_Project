clear
clc

%% 

% Assume that tide doesn't vary much within the PICCOLO survey area
% Choose one PICCOLO station as reference point for the CATS Model
% Check against other PICCOLO stations + reference point further away

%% station coordinates

% coordinates for PIC 1 NW corner
long_sim_1 = -(55 + 10.553/60);
lat_sim_1 = -(64 + 30.924/60);

% coordinates for PIC4 NE corner
long_sim_4 = -(54 + 55.573/60);
lat_sim_4 = -(64 + 31.460/60);

% coordinates for simulation PIC 7 SE corner
long_sim_7=-54.9;
lat_sim_7=-64.6;

% coordinates for PIC 10 SW corner
long_sim_10 = -(55 + 12.268/60);
lat_sim_10 = -(64 + 36.968/60);

% 70 km away (on the shelf)
long_sim_70 = -57.5;
lat_sim_70 = -65.1;


%% map

% PIC7
% convert from degrees to km (polar stereographic coord. 
[x_sim_7,y_sim_7]=xy_ll_CATS2008(long_sim_7,lat_sim_7,'F');


% get bathymetry from cats 
[x,y,H]=tmd_get_bathy('Model_CATS2008');
cats_bathy = figure(1)
clf
pcolor(x,y,H)
shading flat
hold on
scatter(x_sim_7,y_sim_7,101,'filled', 'r'); % location of PIC7
hold off 

% zoomed in 

cats_bathy2 = figure(2)
clf
pcolor(x,y,H)
xlim([0 1000]);
ylim([2500,3000]);
shading flat
hold on
scatter(x_sim_7,y_sim_7,101,'filled', 'r')
hold off 

clearvars x* y* H
%%  Depths

% observed bottom depth from event log
obs_bdepth_1 = 364;
obs_bdepth_4 = 385;
obs_bdepth_7 = 358; 
obs_bdepth_10 = 390;

% observed depth:

lat = ncread('SD025_piccolo_site_0002d.nc', 'y'); % decimal degrees
lon = ncread('SD025_piccolo_site_0002d.nc', 'x');
z = - ncread('SD025_piccolo_site_0002d.nc', 'z'); % Topography (m) z(x,y)

% load coordinates for the four PIC stations
lat1 = ncread('ctd_sd025_069_psal.nc', 'latitude')';
lon1 = ncread('ctd_sd025_069_psal.nc', 'longitude')';
lat4 = ncread('ctd_sd025_071_psal.nc', 'latitude')';
lon4 = ncread('ctd_sd025_071_psal.nc', 'longitude')';
lat7 = ncread('ctd_sd025_075_psal.nc', 'latitude')';
lon7 = ncread('ctd_sd025_075_psal.nc', 'longitude')';
lat10 = ncread('ctd_sd025_081_psal.nc', 'latitude')';
lon10 = ncread('ctd_sd025_081_psal.nc', 'longitude')';

% find where depth coord correspond to station coord
ind_lat_1 = find(lat <= max(lat1) & lat >= min(lat1));
%ind_lon_1 = find(lon <= max(lon1) & lon >= min(lon1)); 
%ind_lat_4 = find(lat <= max(lat4) & lat >= min(lat4));
%ind_lon_4 = find(lon <= max(lon4) & lon >= min(lon4));
ind_lat_7 = find(lat <= max(lat7) & lat >= min(lat7));
ind_lon_7 = find(lon <= max(lon7) & lon >= min(lon7));
ind_lat_10 = find(lat <= max(lat10) & lat >= min(lat10));
ind_lon_10 = find(lon <= max(lon10) & lon >= min(lon10));


% there is no coordinate that corresponds exactly to lon1, lat4 & lon4, so
% we find the nearest index
ind_lon_1 = dsearchn(lon, lon1);
ind_lon_1 = ind_lon_1(1);
ind_lat_4 = dsearchn(lat, lat4);
ind_lat_4 = ind_lat_4(1);
ind_lon_4 = dsearchn(lon, lon4);
ind_lon_4 = ind_lon_4(1);


% find correct depth and take the mean

z_1 = z(ind_lon_1, ind_lat_1);
mean_z_1 = mean(z_1, 'all', 'omitnan');

z_4 = z(ind_lon_4, ind_lat_4);
mean_z_4 = mean(z_4, 'all', 'omitnan');

z_7 = z(ind_lon_7, ind_lat_7);
mean_z_7 = mean(z_7, 'all', 'omitnan');

z_10 = z(ind_lon_10, ind_lat_10);
mean_z_10 = mean(z_10, 'all', 'omitnan');

% difference between bathymetry data and event log depth:

diff_1 = mean_z_1 - obs_bdepth_1;
diff_4 = mean_z_4 - obs_bdepth_4;
diff_7 = mean_z_7 - obs_bdepth_7;
diff_10 = mean_z_10 - obs_bdepth_10;

% 1-2 metres different in depths


% use tmd and CATS to extract depth for our reference latitudes
[~,~,model_depth_1,~]=tmd_extract_HC('Model_CATS2008',lat_sim_1,long_sim_1,'z',[]);
[~,~,model_depth_4,~]=tmd_extract_HC('Model_CATS2008',lat_sim_4,long_sim_4,'z',[]); 
[~,~,model_depth_7,~]=tmd_extract_HC('Model_CATS2008',lat_sim_7,long_sim_7,'z',[]);  
[~,~,model_depth_10,~]=tmd_extract_HC('Model_CATS2008',lat_sim_10,long_sim_10,'z',[]); 
[~,~,model_depth_60,~]=tmd_extract_HC('Model_CATS2008',lat_sim_70,long_sim_70,'z',[]); 

% bedmachine bathymetry for PIC 1
[x_1,y_1] = ll2xy(lat_sim_1,long_sim_1,-1); % -1 (SH). lat, lon to x, y ps
bc_bathymetry_1 = interpBedmachineAntarctica(x_1, y_1, 'bed');

% bedmachine bathymetry for PIC 4
[x_4,y_4] = ll2xy(lat_sim_4,long_sim_4,-1); 
bc_bathymetry_4 = interpBedmachineAntarctica(x_4, y_4, 'bed');

% bedmachine bathymetry for PIC 7
[x_7,y_7] = ll2xy(lat_sim_7,long_sim_7,-1); % 
bc_bathymetry_7 = interpBedmachineAntarctica(x_7, y_7, 'bed');

% bedmachine bathymetry for PIC 10
[x_10,y_10] = ll2xy(lat_sim_10,long_sim_10,-1); % 
bc_bathymetry_10 = interpBedmachineAntarctica(x_10, y_10, 'bed');

% bedmachine bathymetry for 60km away
[x_70,y_70] = ll2xy(lat_sim_70,long_sim_70,-1); % 
bc_bathymetry_70 = interpBedmachineAntarctica(x_70, y_70, 'bed');

% depth differences CATS vs Bedmachine
depth_diff_1 = model_depth_1 - abs(bc_bathymetry_1); %Cats shallower
depth_diff_4 = model_depth_4 - abs(bc_bathymetry_4); %Cats deeper (but not by much)
depth_diff_7 = model_depth_7 - abs(bc_bathymetry_7); %Cats shallower
depth_diff_10 = model_depth_1 - abs(bc_bathymetry_10); % Cats shallower

% generally, CATS is shallower than bedmachine. This is confirmed in the
% script bathymetry_final.m 

% depth differences observation vs CATS

obs_vs_cats_1 = mean_z_1 - model_depth_1;
obs_vs_cats_4 = mean_z_4 - model_depth_4;
obs_vs_cats_7 = mean_z_7 - model_depth_7;
obs_vs_cats_10 = mean_z_10 - model_depth_10;

% CATS is shallower, from 11m for PIC1 to 131m for PIC10

% depth differences observation vs Bedmachine

obs_vs_bc_1 = mean_z_1 - abs(bc_bathymetry_1);
obs_vs_bc_4 = mean_z_4 - abs(bc_bathymetry_4);
obs_vs_bc_7 = mean_z_7 - abs(bc_bathymetry_7);
obs_vs_bc_10 = mean_z_10 - abs(bc_bathymetry_10);

% bc slightly shallower (except for PIC7, where it is 4m deeper)

clearvars -except bc_bathymetry_70 lat_sim* long_sim* model_depth_7 mean* bc_bathymetry_7
%% Get modelled u, v

% dates from start of feb (tmd needs serial dates)
modelled_dates = (15:(0.5/24):35)+datenum(2023,2,0);


% extract tidal harmonics for all tidal constituents combined for our
% stations + further away

[cats_u_1, ConList] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_1, long_sim_1,'u', []);
[cats_v_1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_1, long_sim_1,'v', []);
[cats_U_1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_1, long_sim_1,'U', []);
[cats_V_1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_1, long_sim_1,'V', []);

[cats_u_4, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_4, long_sim_4,'u', []);
[cats_v_4, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_4, long_sim_4,'v', []);
[cats_U_4, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_4, long_sim_4,'U', []);
[cats_V_4, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_4, long_sim_4,'V', []);

[cats_u_7, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'u', []);
[cats_v_7, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'v', []);
[cats_U_7, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'U', []);
[cats_V_7, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'V', []);

[cats_u_10, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_10, long_sim_10,'u', []);
[cats_v_10, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_10, long_sim_10,'v', []);
[cats_U_10, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_10, long_sim_10,'U', []);
[cats_V_10, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_10, long_sim_10,'V', []);

[cats_u_70, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_70, long_sim_70,'u', []);
[cats_v_70, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_70, long_sim_70,'v', []);
[cats_U_70, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_70, long_sim_70,'U', []);
[cats_V_70, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_70, long_sim_70,'V', []);



depthU=cats_U_7./(cats_u_7*1e-2); % *1e-2 due to u being in cm/s
depthV=cats_V_7./(cats_v_7*1e-2); 

% check that the model is giving the right values by checking
% against modelled depth

depthU=cats_U_7./(cats_u_7*1e-2); % *1e-2 due to u being in cm/s
depthV=cats_V_7./(cats_v_7*1e-2); 

tolerance = 1e-9;

if max(abs(model_depth_7(:) -depthU(:))) < tolerance
    disp('Model depth is equal depthU')
    if max(abs(model_depth_7(:) -depthV(:))) < tolerance
        disp('Model depth is equal depthV')
        if max(abs(depthU(:)-depthV(:))) < tolerance
           disp('DepthU is equal depthV')
        end
    end
end 

% depths are equal

% correct model transport against observed depth, to get more accurate
% velocities
modelled_u_1=cats_U_1./mean_z_1; 
modelled_v_1=cats_V_7./mean_z_1;
modelled_u_4=cats_U_4./mean_z_4; 
modelled_v_4=cats_V_4./mean_z_4;
modelled_u_7=cats_U_7./mean_z_7; 
modelled_v_7=cats_V_7./mean_z_7;
modelled_u_10=cats_U_10./mean_z_10; 
modelled_v_10=cats_V_10./mean_z_10;

% 60 km away corrected using bedmachine bathymetry 
modelled_u_60=cats_U_70./abs(bc_bathymetry_70); 
modelled_v_60=cats_V_70./abs(bc_bathymetry_70);


%% extract harmonics for different constituents

% from observations, the tide is diurnal 
% extract the main diurnal constituents: k1, o1, p1

% lunar-solar declinational diurnal constituent K1
[cats_u_k1, conList] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'u', [5]);
[cats_v_k1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'v', [5]);
[cats_U_k1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'U', [5]);
[cats_V_k1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'V', [5]);

% principal lunar diurnal O1
[cats_u_o1, conlist] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'u', [6]);
[cats_v_o1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'v', [6]);
[cats_U_o1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'U', [6]);
[cats_V_o1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'V', [6]);


% principal solar diurnal constituent P1
[cats_u_p1, conList] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'u', [7]);
[cats_v_p1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'v', [7]);
[cats_U_p1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'U', [7]);
[cats_V_p1, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'V', [7]);


% from analysing the time-series with the above consituents. the biggest
% diurnal influence comes from o1 and k1. 

% combining o1 and k1:
[cats_u_d, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'u', [5,6]);
[cats_v_d, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'v', [5,6]);
[cats_U_d, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'U', [5,6]);
[cats_V_d, ~] =  tmd_tide_pred('Model_CATS2008', modelled_dates, lat_sim_7, long_sim_7,'V', [5,6]);


% correct against observed depth (only plotting for PIC7):
modelled_u_d=cats_U_d/mean_z_7; 
modelled_v_d=cats_V_d./mean_z_7;

modelled_u_k1=cats_U_k1/mean_z_7; 
modelled_v_k1=cats_V_k1./mean_z_7;


modelled_u_o1=cats_U_o1/mean_z_7; 
modelled_v_o1=cats_V_o1./mean_z_7;

modelled_u_p1=cats_U_p1/mean_z_7; 
modelled_v_p1=cats_V_p1./mean_z_7;

clearvars -except modelled* *7

%% Read in SADCP data

load PIC 

% depth avg. u and v (ie. timeseries) 
u_timemean=mean(PIC.u(5:30,:), 'omitnan');
v_timemean=mean(PIC.v(5:30,:), 'omitnan');

% residuals 
u_rem = PIC.u - u_timemean;
v_rem = PIC.v - v_timemean;

% time variable (numeric)
obs_time= datenum(PIC.t) -datenum(2023,2,0);
modelled_time=modelled_dates-datenum(2023,2,0);

% time variable (datetime)
o_t = PIC.t;
m_t = datetime(datestr(modelled_dates));


% Tierra de O'Higgins, Península, Antarctica, Moon phase
% New moon: 20th Feb (04:05) Spring tide
% First quarter: 27th Feb (05:05) Neap tide
% Full Moon: 7th March (09:40) Spring tide

% find where our modelled time is equal to these dates
 nm = find(m_t == datetime('2023-02-20 04:00'));
 fq = find(m_t == datetime('2023-02-27 05:00'));
 fm = find(m_t == datetime('2023-03-07')); % last date of time-series

%% CATS plot for PIC7

% plot modelled u,v and observed u,v
CATS_PIC7 = figure(3)
clf
subplot(2,1,1)
plot(o_t,double(u_timemean),'b','linewidth',2); % observed
hold on 
plot(m_t,modelled_u_7,'r','linewidth',2); % modelled

% tidal constituents
plot(m_t,modelled_u_o1,'color', [0.4660 0.7740 0.1880],'linewidth',2) 
plot(m_t,modelled_u_p1,'color', [0.3010 0.7450 0.9330],'linewidth',2) 
plot(m_t,modelled_u_k1,'color', [0.9 0.67 0.2] ,'linewidth',2)
plot(m_t, modelled_u_d, 'color', [1 0.4470 0.7410],'linewidth',2 )

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('Depth averaged u (m/s)','color','k','FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('u (m/s)','FontSize', 22);
xlabel('Date','FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
legend('Observed velocity','Modelled corrected velocity', 'O1', 'P1', 'K1','O1+K1', 'Location', 'southeast', 'FontSize', 10);
hold off

subplot(2,1,2)
plot(o_t,double(v_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_v_7,'r','linewidth',2) % modelled 

% tidal constituents
plot(m_t,modelled_v_o1,'color', [0.4660 0.7740 0.1880],'linewidth',2); 
plot(m_t,modelled_v_p1,'color', [0.3010 0.7450 0.9330],'linewidth',2) 
plot(m_t,modelled_v_k1,'color', [0.9 0.67 0.2],'linewidth',2)
plot(m_t, modelled_v_d, 'color', [1 0.4470 0.7410],'linewidth',2 )

% moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)
title('Depth averaged v (m/s)','color','k', 'FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('v m/s','FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
xlabel('Date','FontSize', 22);
legend('Observed velocity','Modelled corrected velocity', 'O1', 'P1', 'K1','O1+K1', 'Location', 'southeast', 'FontSize', 10);
hold off
set(gcf,'color','w');


%% u 4 locations

% plot modelled u for 4 different locations
CATS_box_u = figure(4)
clf
subplot(2,2,1)
plot(o_t,double(u_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_u_7,'r','linewidth',2) % modelled

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('Depth Averaged u (m/s) for PIC 7','color','k','FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('u (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
legend('Observed velocity','Modelled corrected velocity', 'Location', 'southeast','FontSize', 10);
hold off


subplot(2,2,2)
plot(o_t,double(u_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_u_4,'r','linewidth',2) % modelled

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)


title('Depth Averaged u (m/s) for PIC 4','color','k','FontSize',27)
ax = gca;
ylabel('u (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
ax.FontSize = 22;
legend('Observed velocity','Modelled corrected velocity', 'Location', 'southeast','FontSize', 10);
hold off

subplot(2,2,3)

plot(o_t,double(u_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_u_1,'r','linewidth',2) % modelled

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('Depth Averaged u (m/s) for PIC 1','color','k','FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('u (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
legend('Observed velocity','Modelled corrected velocity', 'Location', 'southeast','FontSize', 10);
hold off

subplot(2,2,4)

plot(o_t,double(u_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_u_10,'r','linewidth',2) % modelled


% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('Depth Averaged u (m/s) for PIC 10','color','k','FontSize',27)
ax = gca;
ylabel('u (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
ax.FontSize = 22;
legend('Observed velocity','Modelled corrected velocity', 'Location', 'southeast', 'FontSize', 10);
hold off

set(gcf,'color','w');

%% v 4 locations

%%% plot modelled v for 4 different locations
CATS_box_v = figure(5)
clf
subplot(2,2,1)
plot(o_t,double(v_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_v_7,'r','linewidth',2) % modelled

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('Depth Averaged v (m/s) for PIC 7','color','k','FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('v (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
legend('Observed velocity','Modelled velocity', 'Location', 'southeast', 'FontSize', 10);
hold off


subplot(2,2,2)
plot(o_t,double(v_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_v_4,'r','linewidth',2) % modelled

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)


title('Depth Averaged v (m/s) for PIC 4','color','k','FontSize',27)
ax = gca;
ylabel('v (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
ax.FontSize = 22;
legend('Observed velocity','Modelled velocity', 'Location', 'southeast', 'FontSize', 10);
hold off

subplot(2,2,3)

plot(o_t,double(v_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_v_1,'r','linewidth',2) % modelled

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('Depth Averaged v (m/s) for PIC 1','color','k','FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('v (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
legend('Observed velocity','Modelled velocity', 'Location', 'southeast', 'FontSize', 10);
hold off

subplot(2,2,4)

plot(o_t,double(v_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_v_10,'r','linewidth',2) % modelled


% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)


title('Depth Averaged v (m/s) for PIC 10','color','k','FontSize',27)
ax = gca;
ylabel('v (m/s)', 'FontSize', 22);
xlabel('Date', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
ax.FontSize = 22;
legend('Observed velocity','Modelled velocity', 'Location', 'southeast', 'Fontsize', 10);
hold off

set(gcf,'color','w');

%% PIC7 vs 70km away 


% plot modelled u,v and observed u,v
CATS_PIC7 = figure(6)
clf
subplot(2,2,1)
plot(o_t,double(u_timemean),'b','linewidth',2); % observed
hold on 
plot(m_t,modelled_u_7,'r','linewidth',2); % modelled

% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('PIC7 u (m/s)','color','k','FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('u (m/s)', 'FontSize', 22);
xlabel('Days', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
legend('Observed velocity','Modelled velocity', 'Location', 'southeast', 'FontSize', 10);
hold off

subplot(2,2,3)
plot(o_t,double(v_timemean),'b','linewidth',2) % observed
hold on 
plot(m_t,modelled_v_7,'r','linewidth',2) % modelled 

% moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('PIC7 v (m/s)','color','k', 'FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('v m/s', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
xlabel('Days', 'FontSize', 22);
legend('Observed velocity','Modelled velocity', 'Location', 'southeast', 'FontSize', 10);
hold off

subplot(2,2,2)
plot(m_t,modelled_u_60,'r','linewidth',2); % modelled


% Moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)

title('70 km away u (m/s)','color','k','FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('u (m/s)', 'FontSize', 22);
xlabel('Days', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
legend('Modelled velocity', 'Location', 'southeast', 'FontSize', 10);

subplot(2,2,4)
plot(m_t,modelled_v_60,'r','linewidth',2) % modelled 

% moon phases 
xline(m_t(nm), 'k', 'New Moon', 'FontSize', 14)
xline(m_t(fq), 'k', 'First Quarter', 'FontSize', 14)
xline(m_t(fm), 'k', 'Full Moon', 'FontSize', 14)
title('70 km away v (m/s)','color','k', 'FontSize',27)
ax = gca;
ax.FontSize = 22;
ylabel('v m/s', 'FontSize', 22);
yticks([-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ylim([-0.5 0.5]);
xlabel('Days', 'FontSize', 22);
legend('Modelled velocity', 'Location', 'southeast', 'FontSize', 10);

set(gcf,'color','w');

%% save modelled values in structure

% use CATS data only in the same time period as PIC.

ind = find(modelled_time >= obs_time(1) & modelled_time <= obs_time(end));

cats_u = modelled_u_7(ind(1):ind(end),1);
cats_v = modelled_v_7(ind(1):ind(end),1);
time = modelled_time(1,ind(1):ind(end));

CATS.u = modelled_u_7;
CATS.u_ot = cats_u;
CATS.v = modelled_v_7;
CATS.v_ot = cats_v;
CATS.lat = lat_sim_7;
CATS.lon = lat_sim_7;
CATS.depth = model_depth_7;
CATS.bc_depth = bc_bathymetry_7;
CATS.obs_depth = mean_z_7;
CATS.dates = m_t;
CATS.datenr = modelled_time;
CATS.time = time;
CATS.u_d = modelled_u_d;
CATS.v_d = modelled_v_d;

save('CATS.mat', 'CATS');
