clc
clear

%% bedmachine

% load in netcdf from NSIDC
bc_H = ncread('BedMachineAntarcticav3.nc', 'bed');

% specify area and resolution 
lat_range = [(-70:0.005:-55)];% lat, lon needs to be same size for ll2xy
lon_range = [(-63:0.005:-48)]; 

% create lat lon meshgrid for interpolation 
[meshlat_range, meshlon_range] = meshgrid(lat_range,lon_range);

% convert meshgrid to polar stereographic coord (lat,lon) = (x,y)
% interpolate bedmachine to meshgrid.
% note: lat, lon needs to be same size
[x_bc,y_bc] = ll2xy(meshlat_range,meshlon_range,-1); % -1 = SH
bc_bathymetry = interpBedmachineAntarctica(x_bc, y_bc, 'bed'); 

clear x_bc y_bc

%% plot bedmachine bathymetry using m_map

% pic mooring location
my = -(64 + 34.380/60); % lat
mx = -(55 + 03.544/60); % lon

b_bc = figure(1)
clf
ax = gca;
m_proj('lambert','lat',[-68 -63],'long',[-62 -48]);
%m_proj('equidistant cylindrical','lat',[-68 -63],'long',[-62 -48]);
m_pcolor(meshlon_range,meshlat_range,bc_bathymetry);
%caxis(ax,[-4500 0])
colormap(ax, 'gray')
m_gshhs_i('patch','r'); %coastline
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
m_northarrow(-61,-63.5, 0.6);
ax1=m_contfbar(0.93,[.5 .9],[-4500 0],[-4500:100:000], 'endpiece','no', 'LineStyle', 'none');
%m_shadedrelief(lon_range, lat_range, bc_bathymetry') % add shading
hold on 
m_scatter(mx,my,101,'filled', 'b', ''); % plot mooring location 
m_text(mx,(my - 0.1),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
ax.FontSize = 22;
ax1.FontSize = 22;
title('Bedmachine Bathymetry', 'FontSize', 27);
xlabel(ax1,'Depth (Meters)','color','k', 'FontSize',22);
xlabel('Longitude', 'FontSize',22);
ylabel('Latitude', 'FontSize', 22);
set(gcf,'color','w');

%% quiver plot of u,v

% load in data from PICCOLO
load PIC

b_bc_quiver = figure(2)
clf
ax = gca;
m_proj('lambert','lat',[-66 -64],'long',[-57 -54]);
%m_proj('equidistant cylindrical','lat',[-65 -64],'long',[-56 -54]);
m_pcolor(meshlon_range,meshlat_range,bc_bathymetry);
caxis(ax,[-4500 0])
colormap(ax,gray)
m_gshhs_i('patch','r'); %coastline
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
m_northarrow(-61,-63.5, 0.6);
ax1=m_contfbar(0.93,[.5 .9],[-4500 0],[-4500:100:000], 'endpiece','no');
m_shadedrelief(lon_range, lat_range, bc_bathymetry') % add shading 
title('Bedmachine Bathymetry (M-Map)', 'FontSize', 27);
xlabel(ax1,'Depth (Meters)','color','k', 'FontSize',22);
xlabel('Longitude', 'FontSize',22);
ylabel('Latitude', 'FontSize', 22);
hold on 
m_quiver(PIC.lon,PIC.lat,mean(PIC.u, 'omitnan'),mean(PIC.v, 'omitnan'),0);
hold off
% 0 removes the grid scaling and scales by the magnitude of u and v
set(gcf,'color','w');
%% plot bedmachine bathymetry using Antarctic Mapping Tools

bc_bath = figure(3)
clf
pcolor(meshlon_range, meshlat_range, bc_bathymetry);
ax = gca;
shading flat;
%colormap(ax, flipud(cbrewer2('BrBG',256)));
cmocean topo % sets colormap to topography
cb = colorbar;
clim([-2000 2000]);
shadem(2,[218 81]), % hillshade 
ylabel(cb,'Elevation above mean sea level (m)', 'FontSize',22)
xlim([-63 -50]);
ylim([-69 -62]);
set(ax,'XTick',[-62 -60 -58 -56 -54 -52 -50],'XTickLabel',{['62', ...
    char(176),'W'],['60',char(176),'W'],['58',char(176),'W'], ...
    ['56',char(176),'W'], ['54',char(176),'W'], ['52',char(176),'W'], ...
    ['50',char(176),'W']});
set(ax,'YTick',[-68 -66 -64 -62],'YTickLabel',{['68',char(176),'S'], ...
    ['66',char(176),'S'],['64',char(176),'S'], ['62',char(176),'S']});
ax.FontSize = 22;
hold on 
scatter(mx,my,101,'filled', 'r', ''); % plot mooring location 
text(mx,(my - 0.1),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
title('Bedmachine Bathymetry', 'FontSize',27);
ylabel('Latitude S', 'FontSize',22);
xlabel('Longitude W', 'FontSize',22);
set(gcf,'color','w');

%% cats 


% extract model bathymetry for our area
[~,~,c_Z,~]=tmd_extract_HC('Model_CATS2008',meshlat_range,meshlon_range, ...
    'z',[]);

% c_Z is postive downwards. Make it negative downwards, similar to bc &
% gebco
cats_bathymetry = -c_Z;


figure(4)
clf
pcolor(meshlon_range, meshlat_range, cats_bathymetry);
ax = gca;
shading flat;
%colormap(ax, flipud(cbrewer2('BrBG',256)));
cmocean topo % sets colormap to topography
cb = colorbar;
clim([-2000 2000]);
shadem(2,[218 81]), % hillshade 
ylabel(cb,'Elevation above mean sea level (m)', 'FontSize',22)
xlim([-63 -50]);
ylim([-69 -62])
set(ax,'XTick',[-62 -60 -58 -56 -54 -52 -50],'XTickLabel',{['62', ...
    char(176),'W'],['60',char(176),'W'],['58',char(176),'W'], ...
    ['56',char(176),'W'], ['54',char(176),'W'], ['52',char(176),'W'], ...
    ['50',char(176),'W']});
set(ax,'YTick',[-68 -66 -64 -62],'YTickLabel',{['68',char(176),'S'], ...
    ['66',char(176),'S'],['64',char(176),'S'], ['62',char(176),'S']});
ax.FontSize = 22;
hold on 
scatter(mx,my,101,'filled', 'r', ''); % plot mooring location 
text(mx,(my - 0.1),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
title('CATS Bathymetry', 'FontSize',27);
ylabel('Latitude S', 'FontSize',22);
xlabel('Longitude W', 'FontSize', 22);
set(gcf,'color','w');

%% cats with M_Map


b_cats = figure(5)
clf
ax = gca;
m_proj('lambert','lat',[-68 -63],'long',[-62 -48]);
%m_proj('equidistant cylindrical','lat',[-68 -63],'long',[-62 -48]);
m_pcolor(meshlon_range,meshlat_range,cats_bathymetry);
caxis(ax,[-4500 0])
colormap(ax, 'gray')
m_gshhs_i('patch','r'); %coastline
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
m_northarrow(-61,-63.5, 0.6);
ax1=m_contfbar(0.93,[.5 .9],[-4500 0],[-4500:100:000], 'endpiece','no', 'LineStyle', 'none');
hold on 
m_scatter(mx,my,101,'filled', 'b', '');
m_text(mx,(my - 0.1),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
ax.FontSize = 22;
ax1.FontSize = 22;
title('CATS Bathymetry', 'FontSize', 27);
xlabel(ax1,'Depth (Meters)','color','k', 'FontSize',22);
xlabel('Longitude', 'FontSize',22);
ylabel('Latitude', 'FontSize', 22);
set(gcf,'color','w');
%% gebco

%The information for ice-surface elevation and under-ice 
% topography/bathymetry is taken from data based on MEaSUREs 
% BedMachine Antarctica, Version 2 (Morlighem, M. et al 2020).


% South of 60°S, the land/ ice-surface elevation topography is largely 
% determined from MEaSUREs BedMachine Antarctica, Version 2 (Morlighem, 
% M. et al 2020). For areas north of 60°N, land data are largely taken 
% from the Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010) 
% data set (Danielson, J.J., and Gesch, D.B., 2011).

% load in variables from the gebco netcdf
geb_lat = ncread('gebco_2023_2.nc', 'lat');
geb_lon = ncread('gebco_2023_2.nc', 'lon');
geb_bathy = double(ncread('gebco_2023_2.nc', 'elevation')); % elevation 
% above mean sea-level

% find where gebco latitude is within specified range
indlat = find(geb_lat > lat_range(1) & geb_lat < lat_range(1,end));
indlon = find(geb_lon > lon_range(1) & geb_lon < lon_range(1,end));

% filter for index
geb_bathy = geb_bathy(:,indlat);

% create gebco meshgrid
[meshlat_geb, meshlon_geb] = meshgrid(geb_lat, geb_lon);

%% gebco m_map

figure(6)
clf
ax = gca;
m_proj('lambert','lat',[-68 -63],'long',[-62 -48]);
%m_proj('equidistant cylindrical','lat',[-68 -63],'long',[-62 -48]);
m_pcolor(meshlon_geb,meshlat_geb,geb_bathy);
caxis(ax,[-4500 0])
colormap(ax,gray)
m_gshhs_i('patch','r'); %coastline
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
m_northarrow(-61,-63.5, 0.6);
ax1=m_contfbar(0.93,[.5 .9],[-450000 0],[-4500:100:000], 'endpiece','no', ...
    'LineStyle', 'none');
ax.FontSize = 22;
ax1.FontSize = 22;
hold on 
m_scatter(mx,my,101,'filled', 'b', ''); % plot mooring location 
m_text(mx,(my - 0.1),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
title('GEBCO Bathymetry', 'FontSize', 27);
xlabel(ax1,'Depth (Meters)','color','k', 'FontSize',22);
xlabel('Longitude', 'FontSize',22);
ylabel('Latitude', 'FontSize', 22);
set(gcf,'color','w');

%% gebco bathymetry using Antarctic Mapping Tools

figure(7)
clf
pcolor(meshlon_geb, meshlat_geb, geb_bathy);
ax = gca;
shading flat;
%colormap(ax, flipud(cbrewer2('BrBG',256)));
cmocean topo % sets colormap to topography
cb = colorbar;
clim([-2000 2000]);
shadem(2,[218 81]), % hillshade 
ylabel(cb,'Elevation above mean sea level (m)', 'FontSize',22)
xlim([-63 -50]);
ylim([-69 -62])
set(ax,'XTick',[-62 -60 -58 -56 -54 -52 -50],'XTickLabel',{['62', ...
    char(176),'W'],['60',char(176),'W'],['58',char(176),'W'], ...
    ['56',char(176),'W'], ['54',char(176),'W'], ['52',char(176),'W'], ...
    ['50',char(176),'W']});
set(ax,'YTick',[-68 -66 -64 -62],'YTickLabel',{['68',char(176),'S'], ...
    ['66',char(176),'S'],['64',char(176),'S'], ['62',char(176),'S']});
ax.FontSize = 22;
hold on 
scatter(mx,my,101,'filled', 'r', ''); % plot mooring location 
text(mx,(my - 0.1),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
title('GEBCO Bathymetry', 'FontSize',27);
ylabel('Latitude S', 'FontSize',22);
xlabel('Longitude W', 'FontSize',22);
set(gcf,'color','w');
%% difference maps 

% the difference between bedmachine and CATS
bc_diff_cats = bc_bathymetry - cats_bathymetry;
% negative = bc deeper

% gebco has not got the same resolution as bedmachine or cats, so we need 
% to interpolate on top of bedmachine grid using linear interpolation:
geb_bathymetry_interp = interp2(geb_lon', geb_lat', geb_bathy, ...
    meshlon_range, meshlat_range);

% difference between gebco and bedmachine
geb_diff_bc = geb_bathymetry_interp' - bc_bathymetry;
% negative = gebco deeper


% difference between gebco and cats
geb_diff_cats = geb_bathymetry_interp' - cats_bathymetry;
% negative = gebco deeper

figure(8)
clf
subplot(1,3,1)
pcolor(meshlon_range, meshlat_range, bc_diff_cats)
ax = gca;
title('Bedmachine - CATS', 'FontSize',27)
shading interp;
caxis([-300 300])
cmocean balance
cb = colorbar;
ylabel(cb,'Difference (m)', 'FontSize',22)
xlim([-63 -54]);
ylim([-69 -62])
set(ax,'XTick',[-62 -60 -58 -56 -54],'XTickLabel',{['62', ...
    char(176),'W'],['60',char(176),'W'],['58',char(176),'W'], ...
    ['56',char(176),'W'], ['54',char(176),'W']});
set(ax,'YTick',[-68 -66 -64 -62],'YTickLabel',{['68',char(176),'S'], ...
    ['66',char(176),'S'],['64',char(176),'S'], ['62',char(176),'S']});
ax.FontSize = 22;
%ylabel('Latitude ', 'FontSize',22)
%xlabel('Longitude', 'FontSize',22);

subplot(1,3,2)
pcolor(meshlon_range, meshlat_range, geb_diff_bc)
ax = gca;
title('GEBCO - Bedmachine', 'FontSize',27)
shading interp;
caxis([-300 300])
cmocean balance
cb = colorbar;
ylabel(cb,'Difference (m)', 'FontSize',22);
xlim([-63 -54]);
ylim([-69 -62])
set(ax,'XTick',[-62 -60 -58 -56 -54],'XTickLabel',{['62', ...
    char(176),'W'],['60',char(176),'W'],['58',char(176),'W'], ...
    ['56',char(176),'W'], ['54',char(176),'W']});
set(ax,'YTick',[-68 -66 -64 -62],'YTickLabel',{['68',char(176),'S'], ...
    ['66',char(176),'S'],['64',char(176),'S'], ['62',char(176),'S']});
ax.FontSize = 22;
%ylabel('Latitude', 'FontSize',22)
%xlabel('Longitude', 'FontSize',22);

subplot(1,3,3)
pcolor(meshlon_range, meshlat_range, geb_diff_cats)
ax = gca;
title('GEBCO - CATS', 'FontSize',27)
shading interp;
caxis([-300 300])
cmocean balance
cb = colorbar;
ylabel(cb,'Difference (m)', 'FontSize',22)
xlim([-63 -54]);
ylim([-69 -62]);
set(ax,'XTick',[-62 -60 -58 -56 -54],'XTickLabel',{['62', ...
    char(176),'W'],['60',char(176),'W'],['58',char(176),'W'], ...
    ['56',char(176),'W'], ['54',char(176),'W']});
set(ax,'YTick',[-68 -66 -64 -62],'YTickLabel',{['68',char(176),'S'], ...
    ['66',char(176),'S'],['64',char(176),'S'], ['62',char(176),'S']});
ax.FontSize = 22;
%ylabel('Latitude', 'FontSize',22);
%xlabel('Longitude', 'FontSize',22);

set(gcf,'color','w');

%% observed bathymetry

% observed depth:
obs_lat = ncread('SD025_piccolo_site_0002d.nc', 'y'); % decimal degrees
obs_lon = ncread('SD025_piccolo_site_0002d.nc', 'x');
z = ncread('SD025_piccolo_site_0002d.nc', 'z'); % Topography (m) z(x,y)

% find coordinate end-points
lat_range_obs = [max(obs_lat)  min(obs_lat)];
lon_range_obs = [max(obs_lon)  min(obs_lon)];

% find where gebco is within this range
indlat2 = find(geb_lat >= lat_range_obs(end) & geb_lat <= ...
    lat_range_obs(1));
indlon2 = find(geb_lon >= lon_range_obs(end) & geb_lon <= ...
    lon_range_obs(1));

% get the correct lat, lon & bathymetry
gebco_lat_pic = geb_lat(indlat2,1);
gebco_lon_pic = geb_lon(indlon2,1);
geb_bath_pic = geb_bathy(indlon2,indlat2);

% find where bedmachine is whitin this range
indlat3 = find(lat_range >= lat_range_obs(end) & lat_range <= ...
    lat_range_obs(1));
indlon3 = find(lon_range >= lon_range_obs(end) & lon_range <= ...
    lon_range_obs(1));

bc_lat_pic = lat_range(1,indlat3);
bc_lon_pic = lon_range(1,indlon3);
bc_bath_pic = bc_bathymetry(indlon3,indlat3);

% find where CATS is within this range 
% convert the chosen and observed range to polar stereographic coordinates
[xcat,ycat]=xy_ll_CATS2008(lat_range,lon_range, 'F');
[xobs,yobs]=xy_ll_CATS2008(lat_range_obs,lon_range_obs, 'F');  

indlat4 = find(xcat >= xobs(end) & xcat <= xobs(1));
indlon4 = find(ycat >= yobs(end) & ycat <= yobs(1));


cats_lat_pic = lat_range(1,indlat4);
cats_lon_pic = lon_range(1,indlon4);
cats_bath_pic= c_Z(indlon4,indlat4);

%[~,~,c_Z_PIC,~]=tmd_extract_HC('Model_CATS2008',obs_lat,obs_lon,'z',[]); % array sizes incomptatible

%% plot of observations

[meshlon_obs, meshlat_obs] = meshgrid(obs_lon, obs_lat);

figure(9)
clf
ax = gca;
%m_proj('lambert','lat',[-(64 + 48/60) -(64 + 12/60)],'long',[-(55 + 30/60)
%  -(54 + 30/60)]);
m_proj('equidistant cylindrical','lat',[-(64 + 48/60) -(64 + 18/60)], ...
    'long',[-(55 + 30/60) -(54 + 30/60)]);
m_pcolor(meshlon_obs,meshlat_obs,z');
colormap(ax, flipud(cbrewer2('Spectral', 256)));
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
ax1=m_contfbar(0.93,[.5 .9],[-500 0],[-450:1:-300], 'endpiece','no', ...
    'LineStyle','none');
hold on 
m_scatter(mx,my,101,'filled', 'b', ''); % plot mooring location 
m_text(mx,(my-0.01),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
ax.FontSize = 22;
ax1.FontSize = 22; % contourbar 
title('Observed Bathymetry', 'FontSize', 27);
xlabel(ax1,'Depth below mean sealevel (m)','color','k', 'FontSize',22);
set(gcf,'color','w');

%% diff between obs and models

% interpolate bedmachine onto obs grid
bc_obs_interp = interp2(bc_lon_pic, bc_lat_pic, bc_bath_pic', ...
    meshlon_obs, meshlat_obs);

% interpolate cats onto obs grid
cats_obs_interp = interp2(lon_range, lat_range, cats_bathymetry', ...
    meshlon_obs, meshlat_obs);

% interpolate gebco onto obs grid
geb_obs_interp = interp2(gebco_lon_pic', gebco_lat_pic', geb_bath_pic', ...
    meshlon_obs, meshlat_obs);

%% differences from observations

% cats and obs
cats_diff_obs = cats_obs_interp' - z;

% bedmachine and obs
bc_diff_obs = bc_obs_interp' - z;

% gebco and obs 
geb_diff_obs = geb_obs_interp' -z;


%% plot of observations and differences in models 


figure(10)
clf
subplot(2,2,1)
ax = gca;
m_proj('lambert','lat',[-(64 + 48/60) -(64 + 18/60)],'long',[-(55 + 30/60) ...
    -(54 + 30/60)]);
%m_proj('equidistant cylindrical','lat',[-(64 + 48/60) -(64 + 18/60)], ...
%    'long',[-(55 + 30/60) -(54 + 30/60)]);
m_pcolor(meshlon_obs,meshlat_obs,z');
ax1 = colorbar;
ax1.FontSize = 22;
colormap(ax, flipud(cbrewer2('Spectral', 256)));
caxis([-450 -300]);
m_grid('xtick',7,'tickdir','out','fontsize',16, 'box', ('fancy'));
hold on 
m_scatter(mx,my,101,'filled', 'k', ''); % plot mooring location 
m_text(mx,(my - 0.01),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
ax.FontSize = 22;
title('Observations', 'FontSize', 27);
xlabel(ax1,'Depth below mean sea-level (m)','color','k', 'FontSize',22);

subplot(2,2,2)
ax = gca;
m_proj('lambert','lat',[-(64 + 48/60) -(64 + 18/60)],'long',[-(55 + 30/60)
  -(54 + 30/60)]);
%m_proj('equidistant cylindrical','lat',[-(64 + 48/60) -(64 + 18/60)], ...
%    'long',[-(55 + 30/60) -(54 + 30/60)]);
m_pcolor(meshlon_obs,meshlat_obs, cats_diff_obs');
ax1 = colorbar;
colormap(ax, cbrewer2('RdBu', 256));
caxis([-100 100]);
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
hold on 
m_scatter(mx,my,101,'filled', 'k', ''); % plot mooring location 
m_text(mx,(my - 0.01),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
ax.FontSize = 22;
ax1.FontSize = 22; % contourbar 
title('CATS - Observations', 'FontSize', 27);
xlabel(ax1,'Difference (m)','color','k', 'FontSize',22);

subplot(2,2,3)
ax = gca;
m_proj('lambert','lat',[-(64 + 48/60) -(64 + 18/60)],'long',[-(55 + 30/60)
  -(54 + 30/60)]);
%m_proj('equidistant cylindrical','lat',[-(64 + 48/60) -(64 + 18/60)], ...
%    'long',[-(55 + 30/60) -(54 + 30/60)]);
m_pcolor(meshlon_obs,meshlat_obs, bc_diff_obs');
ax1 = colorbar;
colormap(ax, cbrewer2('RdBu', 256));
caxis([-100 100]);
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
hold on 
m_scatter(mx,my,101,'filled', 'k', ''); % plot mooring location 
m_text(mx,(my - 0.01),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
ax.FontSize = 22;
ax1.FontSize = 22; % contourbar 
title('Bedmachine - Observations', 'FontSize', 27);
xlabel(ax1,'Difference (m)','color','k', 'FontSize',22);

subplot(2,2,4)
ax = gca;
m_proj('lambert','lat',[-(64 + 48/60) -(64 + 18/60)],'long',[-(55 + 30/60)
  -(54 + 30/60)]);
%m_proj('equidistant cylindrical','lat',[-(64 + 48/60) -(64 + 18/60)], ...
%    'long',[-(55 + 30/60) -(54 + 30/60)]);
m_pcolor(meshlon_obs,meshlat_obs, geb_diff_obs');
ax1 = colorbar;
colormap(ax, cbrewer2('RdBu', 256));
caxis([-100 100]);
m_grid('xtick',5,'tickdir','out','fontsize',16, 'box', ('fancy'));
hold on 
m_scatter(mx,my,101,'filled', 'k', ''); % plot mooring location 
m_text(mx,(my - 0.01),'Mooring', 'FontSize', 12, 'vertical', 'top');
hold off
ax.FontSize = 22;
ax1.FontSize = 22; % contourbar 
title('GEBCO - Observations', 'FontSize', 27);
xlabel(ax1,'Difference (m)','color','k', 'FontSize',22);

set(gcf,'color','w');
