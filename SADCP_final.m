clear
clc

%% SADCP DATA: 75kHz & 150kHz

% The 75kHz and the 150kHz SADCP did not run simultanously. Due to the bad
% quality of the 75kHz data, it was only run for the beginning of the
% cruise and tested once towards the end. The last SADCP to run was the
% 150kHz.


% Trial cruise timeline from event-log (based on CTD casts):
% Start: 01/02-23 13:18
% End: 09/03-23 19:37 


% PIC timeline from event-log (based on CTD casts):
% Start: Event 87, Cast 66. Coord: 65 12.982 S,  56 01.602 W. Date: 23/02/15 13:52
% YOYO1: Event 108, Cast 85-88. Coord: 64 34.015 S, 55 00.704W. Date: 23/02/18 01:55-03:03	 
% YOYO2: Event 109, Cast 89-156. Coord: 64 34.636 S,  54 58.513 W. Date: 23/02/18 04:19 -  23/02/19 0146 (Moves in a southward direction, but back and forth between E and W)	 
% End: Event 117, Cast 160. Coord: 64 31.414, S 54 55.757 W. Date: 23/02/19 18:57

% PICCOLO data will be referred to as PIC data
%% time

% read in netcdfs
d_days_75 = ncread('os75.nc', 'time'); % units: decimal day
d_days_150 = ncread('os150.nc', 'time');


% convert to duration array
d_days_75 = days(d_days_75);
d_days_150 = days(d_days_150);

% days since 2023-01-01 00:00:00. 
time_75 = datetime('2023-01-01 00:00:00') + d_days_75; 
time_150 = datetime('2023-01-01 00:00:00') + d_days_150;


% check that the time is correct:
disp(time_75(1)); % start: 1-Feb-2023 08:51:47 (before first CTD cast)
disp(time_150(end)); % end: 17-Mar-2023 18:49:43 (return after CTD casts)



%% u,v, depth, lat, lon 

% read in netcdfs
u_75 = ncread('os75.nc', 'u'); 
v_75 = ncread('os75.nc', 'v'); 
depth_75 = ncread('os75.nc', 'depth'); 
u_150 = ncread('os150.nc', 'u');
v_150= ncread('os150.nc', 'v'); 
depth_150 = ncread('os150.nc', 'depth');  

% missing value = 9.999999680285692e+37
nan_val = 9.999999680285692e+37;

% use written function to replace value with NaN
u_75 = replaceValueWithNaN(u_75, nan_val);
v_75 = replaceValueWithNaN(v_75, nan_val);
depth_75 = replaceValueWithNaN(depth_75, nan_val);
u_150 = replaceValueWithNaN(u_150, nan_val);
v_150 = replaceValueWithNaN(v_150, nan_val);
depth_150 = replaceValueWithNaN(depth_150, nan_val);

clearvars nan_val

% read in netcdfs
lon_75 = ncread('os75.nc', 'lon');
lat_75 = ncread('os75.nc', 'lat'); 
lon_150 = ncread('os150.nc', 'lon'); 
lat_150 = ncread('os150.nc', 'lat'); 


% missing value = 1e+38
nan_val = 1e+38;

lat_75 = replaceValueWithNaN(lat_75, nan_val);
lon_75 = replaceValueWithNaN(lon_75, nan_val);
lat_150 = replaceValueWithNaN(lat_150, nan_val);
lon_150 = replaceValueWithNaN(lon_150, nan_val);


% all data flags are removed by replacing missing values

clear hoursOfDay nan_val nan_val i 

%% save Sir David Attenborough trial cruise data

sda.time75 = time_75;
sda.time150 = time_150;
sda.lat75 = lat_75;
sda.lat150 = lat_150;
sda.lon75 = lon_75;
sda.lon150 = lon_150;
sda.v75 = v_75;
sda.v150 = v_150;
sda.u75 = u_75;
sda.u150 = u_150;
sda.depth75 = depth_75;
sda.depth150 = depth_150;


save('sda.mat', 'sda');
%% Pcolor plot for current velocities the whole cruise

% transpose time to row vectors
t_time_75 = time_75';
t_time_150 = time_150';


% period matrices

wholeperiod150 = repmat(t_time_150(1, :), 55, 1); % repmat repeats the first entry of each column 55 times in one set of rows for each 
wholeperiod75 = repmat(t_time_75(1,:),55,1);

cruise = figure(1)
clf
ax1 = subplot(2,1,1)
pcolor(t_time_150, depth_150, v_150) 
hold on 
pcolor(t_time_75, depth_75, v_75)
shading interp % smooth
ax1.Colormap = cbrewer2('RdYlBu', 256);
ax1.FontSize = 22;
c = colorbar;
c.Ticks = [-0.5 -0.25 0 0.25 0.5];
c.FontSize = 22;
c.Label.String = 'v (m/s)';
caxis([-1 1]*0.5);
title(' Trial Cruise SADCP v (m/s)', 'FontSize', 27);
xlabel('Date','FontSize',20);
ylabel('Depth (m)', 'FontSize', 22);
xlim([t_time_75(1), t_time_150(end)])
ylim([0, 700])
set(gca, 'YDir','reverse')
hold off


ax2 = subplot(2,1,2)
pcolor(t_time_150, depth_150, u_150);
hold on 
pcolor(t_time_75, depth_75, u_75);
shading interp % smooth
ax2.Colormap = cbrewer2('RdYlBu', 256);
ax2.FontSize = 22;
c = colorbar;
c.Ticks = [-0.5 -0.25 0 0.25 0.5];
c.FontSize = 22;
c.Label.String = 'u (m/s)';
caxis([-1 1]*0.5);
title('Trial Cruise SADCP u (m/s)', 'FontSize',27);
xlabel('Date', 'FontSize',22);
ylabel('Depth (m)','FontSize', 20);
xlim([t_time_150(1), t_time_150(end)]);
ylim([0, 700]);
set(gca, 'YDir','reverse'); % reverse y-axis
hold off
set(gcf,'color','w');

clear sda
%% filter out PICCOLO cruise timeline

PIC_start = datetime('2023-02-15 13:52:00');
PIC_end = datetime('2023-02-19 18:57:00');


% 75kHz
index_75 = find(time_75 >= PIC_start & time_75 <= PIC_end);
PIC_period_75= time_75(index_75);

% 150kHz
index_150 = find(time_150 >= PIC_start & time_150 <= PIC_end);
PIC_period_150 = time_150(index_150);


%% filter u,v,depth, lat, lon 

d_PIC_days75 = d_days_75(index_75(1):index_75(end,1)); % decimal days
u_PIC_75 = u_75(1:55,index_75(1):index_75(end)); 
v_PIC_75 = v_75(1:55,index_75(1):index_75(end));
depth_PIC_75 = depth_75(1:55,index_75(1):index_75(end));
lat_PIC_75 = lat_75(index_75(1):index_75(end),1);
lon_PIC_75= lon_75(index_75(1):index_75(end),1);

d_PIC_days150 = d_days_150(index_150(1):index_150(end,1)); % decimal days
u_PIC_150 = u_150(1:55,index_150(1):index_150(end));
v_PIC_150 = v_150(1:55,index_150(1):index_150(end));
depth_PIC_150 = depth_150(1:55,index_150(1):index_150(end));
lat_PIC_150 = lat_150(index_150(1):index_150(end),1);
lon_PIC_150 = lon_150(index_150(1):index_150(end),1);

%% save PICCOLO data in structure

% 75kHz and 150kHz not combined 

PICCOLO.time75 = PIC_period_75;
PICCOLO.time150 = PIC_period_150;
PICCOLO.lat75 = lat_PIC_75;
PICCOLO.lat150 = lat_PIC_150;
PICCOLO.lon75 = lon_PIC_75;
PICCOLO.lon150 = lon_PIC_150;
PICCOLO.v75 = v_PIC_75;
PICCOLO.v150 = v_PIC_150;
PICCOLO.u75 = u_PIC_75;
PICCOLO.u150 = u_PIC_150;
PICCOLO.depth75 = depth_PIC_75;
PICCOLO.depth150 = depth_PIC_150;
PICCOLO.decimaldays75 = d_PIC_days75;
PICCOLO.decimaldays150 = d_PIC_days150;

save('PICCOLO.mat', 'PICCOLO');

%% Pcolor plot of current velocities for PICCOlO

% transpose time to row vectors
PIC_per_150 = PIC_period_150';
PIC_per_75 = PIC_period_75';


piccolo = figure(2)
clf
ax1 = subplot(2,1,1)
pcolor(PIC_per_150, depth_PIC_150, v_PIC_150) % 150kHz
hold on 
pcolor(PIC_per_75, depth_PIC_75, v_PIC_75) % 75kHz
shading interp % smooth
hold off
ax1.FontSize = 16;
ax1.Colormap = cbrewer2('RdYlBu', 256);
c = colorbar;
c.Ticks = [-0.5 -0.25 0 0.25 0.5];
c.Label.String = 'v (m/s)';
c.FontSize = 22;
caxis([-1 1]*0.5)
set(gca, 'YDir','reverse') % reverse y-axis
xlim([PIC_per_75(1) PIC_per_150(end)]);
ylim([0 400]); % we have no data below 400m
title('PICCOLO SADCP v (m/s)', 'FontSize',27)
xlabel('Date', 'FontSize',22);
ylabel('Depth (m)', 'FontSize',22);

ax2 = subplot(2,1,2)
pcolor(PIC_per_150, depth_PIC_150, u_PIC_150); % 150kHz
hold on 
pcolor(PIC_per_75, depth_PIC_75, u_PIC_75); % 75kHz
shading interp
hold off
ax2.FontSize = 16;
ax2.Colormap = cbrewer2('RdYlBu', 256);
c = colorbar;
c.Ticks = [-0.5 -0.25 0 0.25 0.5];
c.Label.String = 'u (m/s)';
c.FontSize = 22;
caxis([-1 1]*0.5);
set(gca, 'YDir','reverse')
title(' PICCOLO SADCP u');
ylim([0, 400]);
xlim([PIC_per_75(1) PIC_per_150(end)]); % extend x-axis
title('PICCOLO SADCP u (m/s)', 'FontSize',27)
xlabel('Date', 'FontSize',22);
ylabel('Depth (m)', 'FontSize',22);
set(gcf,'color','w');

%% combining 75kHz and 150kHz data for PICCOLO

% The SADCP 75khz and 150khz did not run at the same time.
% We can (with some confidence) combine the two datasets, as we will then 
% have data covering a larger temporal range, and since we are assuming 
% that we have a time-series. This would be beneficial. 

% 75kHz start and end 
start_75 = PIC_period_75(1,1); % 15 Feb 13:56:31
end_75 = PIC_period_75(end,1); % 19 Feb 10:21:00

% 150kHz start and end 
start_150 = PIC_period_150(1,1); % 16 Feb 14:34:10
end_150 = PIC_period_150(end,1); % 19 Feb 18:53:51


% variable matrices
u = [u_PIC_75 u_PIC_150];
v = [v_PIC_75 v_PIC_150];
depth = [depth_PIC_75 depth_PIC_150];
lat = [lat_PIC_75' lat_PIC_150'];
lon = [lon_PIC_75' lon_PIC_150'];
t = [PIC_per_75 PIC_per_150];
dt = [d_PIC_days75' d_PIC_days150']; % decimal day

% sort time in ascending order and get new indices
[time, indx] = sort(t);

% sort variables by the new indices
u = u(:,indx);
v = v(:,indx);
depth = depth(:,indx);
lat = lat(:,indx);
lon = lon(:,indx);
dt = dt(:,indx);

clearvars -except u v depth lat lon time dt

% Now we have all the data we need from the PICCOLO cruise

%% plot PICCOLO data by timeofday

% convert decimal days from duration array to double
d = days(dt);

% floor to retain decimal values corresponding to time of day
td = d - floor(d); 

% sort the data by time of day
[hd,idx] = sort(td);

% arrange variables by the new index
u_sorted = u(:,idx);
v_sorted = v(:,idx);
d_sorted = depth(:,idx);
lat_sorted = lat(:,idx);
lon_sorted = lon(:,idx);
depth_sorted = depth(:,idx);

% create a time matrix of same size as depth, u,v
t_matrix = repmat(hd(1,:),55, 1); % sorted version

% scatter plot by time of day
sctrOneDay = figure(3)
clf
ax1 = subplot(2,1,1)
scatter(t_matrix(:),depth_sorted(:),40,v_sorted(:),'filled'); % (:) vector
ylim([0, 400]); 
hold on 
ax1.Colormap = cbrewer2('RdYlBu', 256);
ax1.FontSize = 22;
c = colorbar;
c.Ticks = [-0.5 -0.25 0 0.25 0.5];
c.Label.String = 'v (m/s)';
c.FontSize = 22;
caxis([-1 1]*0.5);
set(gca, 'YDir','reverse')
hold off
title(' Meridional current velocity', 'FontSize', 27);
xlabel('Time of Day (decimal day)', 'FontSize', 22);
ylabel('Depth (m)', 'FontSize',22);

ax2 = subplot(2,1,2)
scatter(t_matrix(:),depth_sorted(:),40,u_sorted(:),'filled');
ylim([0, 400]); 
hold on 
ax2.Colormap = cbrewer2('RdYlBu', 256);
ax2.FontSize = 22;
c = colorbar;
c.Ticks = [-0.5 -0.25 0 0.25 0.5];
c.FontSize = 22;
c.Label.String = 'u (m/s)';
caxis([-1 1]*0.5);
set(gca, 'YDir','reverse')
hold off
title(' Zonal current velocity', 'FontSize', 27);
xlabel('Time of Day (decimal day)', 'FontSize', 22);
ylabel('Depth (m)', 'FontSize',22);
set(gcf,'color','w');


clearvars ax* d idx c 

%% Scatter plot of mean depth velocity in lat,lon

mean_vel = figure(4)
clf
subplot(1,3,1)
geoscatter(lat,lon,100,mean(v, 'omitnan'), 'filled'); 
grid on
geolimits(gca, [-64.7 -64.35], [-55.3 -55])
ax = gca; % set axis handle
ax.Basemap = 'bluegreen'
ax.FontSize = 22;
c = colorbar;
c.FontSize = 22;
caxis([-0.5, 0.5]);
c.Label.String = 'Mean v (m/s)';
ax.Colormap = cbrewer2('PRGn', 256);
title('Depth avg. meridional velocity', 'FontSize',27);
hold off

subplot(1,3,2)
geoscatter(lat,lon,100,mean(u, 'omitnan'), 'filled'); 
grid on
geolimits(gca, [-64.7 -64.35], [-55.3 -55])
ax = gca; % set axis handle
ax.Basemap = 'bluegreen'
ax.FontSize = 22;
c = colorbar;
c.FontSize = 22;
caxis([-0.5, 0.5]);
c.Label.String = 'Mean u (m/s)';
ax.Colormap = cbrewer2('PRGn', 256);
title('Depth avg. zonal velocity', 'FontSize',27);
hold off

subplot(1,3,3)
geoscatter(lat,lon,100, td, 'filled'); % time of day
grid on
geolimits(gca, [-64.7 -64.35], [-55.3 -55])
ax = gca;
ax.Basemap = 'bluegreen'
ax.FontSize = 22;
c = colorbar;
c.FontSize = 22;
c.Label.String = 'Timeofday';
ax.Colormap = cbrewer2('RdYlGn', 256);
title('Time of day', 'FontSize',27);
hold off
set(gcf,'color','w');

clearvars ax* c gb
%% Save new PIC structure

PIC.t = time;
PIC.decimalday = dt;
PIC.floordd = td;
PIC.u = u;
PIC.v = v;
PIC.depth = depth;
PIC.lat = lat;
PIC.lon = lon;

save('PIC.mat', 'PIC')

%% scatter avg. u, v by time of day


sctvel = figure(5)
clf
scatter(mean(u_sorted, 'omitnan'),mean(v_sorted, 'omitnan'),150, hd, 'filled'); 
grid on 
ax = gca;
ax.FontSize = 22;
title('Depth averaged current velocities', 'FontSize',27);
xlabel('u (m/s)', 'FontSize', 22);
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
ylabel('v (m/s)', 'FontSize', 22);
set(gca,'color','[.5 .7 .9]')
colormap(cbrewer2('RdYlGn', 256));
c = colorbar;
c.FontSize = 2
c.Label.String = 'Time of Day (Decimal Day)';
set(gcf,'color','w');