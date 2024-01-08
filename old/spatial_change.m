
cd('/Users/raphael/Library/CloudStorage/Box-Box/BMM-US/')
load('data/density/inference-trans.mat','g','radar')
addpath(genpath('./functions/'))
cities = readtable('eBird/cities.csv');

%% Read data
q.zone="USA";
q.obs=0;
[res_check, res_obs,q] = queryERD(q);
%save('data/eBird/entireUS.mat','res_check','q')
load('data/eBird/entireUS.mat')

%%
cities.checklist_count(:)=nan;
for i_c=1:height(cities)
    id = res_check.latitude>(cities.latitude(i_c)-.75) &...
        res_check.latitude<(cities.latitude(i_c)+.75)&...
        res_check.longitude>(cities.longitude(i_c)-.75) &...
        res_check.longitude<(cities.longitude(i_c)+.75);

    cities.checklist_count(i_c) = sum(id);
end

thr =cities.checklist_count>.5e5;
figure; hold on;
scatter(cities.longitude(thr),cities.latitude(thr),[],cities.checklist_count(thr),'filled')
text(cities.longitude(thr),cities.latitude(thr),cities.city(thr))
borders('states','k'); colorbar;
axis equal ; axis([min(A.lon) max(A.lon) min(A.lat) max(A.lat)])


%% 

dlatlon=.25;
res_check.lat = round(res_check.latitude/dlatlon)*dlatlon;
res_check.lon = round(res_check.longitude/dlatlon)*dlatlon;

A=groupcounts(res_check,{'lat','lon'});
lon=min(A.lon):dlatlon:max(A.lon);
lat=min(A.lat):dlatlon:max(A.lat);
[LAT,LON]=meshgrid(lat,lon);
M = nan(numel(lat), numel(lon));
A.lon_id = (A.lon-min(A.lon)) / dlatlon+1;
A.lat_id = (A.lat-min(A.lat)) / dlatlon+1;
A.id = sub2ind(size(M),A.lat_id,A.lon_id);

M(round(A.id)) = A.GroupCount; 

figure('position',[0 0 1650 1200]);  set(gcf, 'color', 'none');
tiledlayout(1,1,'TileSpacing','none','Padding','tight')
hold on; xticks([]); yticks([]);
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off; 
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
M(M==0)=nan;
im1=geoshow(LAT,LON,(M'),'DisplayType','texturemap');
colormap(copper)

bordersm('states','color',.3.*[1 1 1]); 
% cityname = ["Kansas City","Chicago","Tompkins","Atlanta","Austin","Boston","Washington","Denver","Minneapolis", "Cincinnati"];
id=ismember(cities.city,cityname);
circlem(cities.latitude(id),cities.longitude(id),111*.75*ones(sum(id),1),'edgecolor',[0.8080    0.1775    0.0166],'linewidth',2)
textm(cities.latitude(id)+1.2,cities.longitude(id),cities.city(id),'HorizontalAlignment','center','color',[0.8080    0.1775    0.0166],'fontsize',18)
circlem(42.4553, -76.4772 ,111*.75,'edgecolor',[0.8080    0.1775    0.0166],'linewidth',2)
textm(42.4553+1.2,-76.4772,"Ithaca",'HorizontalAlignment','center','color',[0.8080    0.1775    0.0166],'fontsize',18)

scatterm(radar.lat, radar.lon,80,'.y')
caxis([0 30000])
colorbar
exportgraphics(gcf, "figures/presentation/ebird_map_5.png",'BackgroundColor','k')

%%

% Combines records of all species. 
% res_obs_check = groupsummary(res_obs(res_obs.species_code=="amerob",:),"sampling_event_id","sum","obs_count");
% res_obs_check = res_obs(res_obs.species_code=="daejun",:);
res_obs_check = res_obs(ismember(res_obs.species_code,speNoc.SPECIES_CODE),:);
 
T = outerjoin(res_check,res_obs_check);
% Ts = groupsummary(T,'species_code',"sum",{'obs_count','effort_hrs'});

% Account for zero: T.obs_count -> nan = not recorded, 0: X (presnce only)
% set x to mean of the group
G=findgroups(T(:, ismember(T.Properties.VariableNames, {'day_of_year','lat','lon','species_code'})));
obs_count_GROUP_Not0 = splitapply(@(x) mean(x(x>0)), T.obs_count, G);
T.obs_count_0 = T.obs_count;
T.obs_count_0(T.obs_count==0) = obs_count_GROUP_Not0(G(T.obs_count==0));
% set absence data to zero for complete list
T.obs_count_0(isnan(T.obs_count)&T.all_obs_reported) = 0;

% Aggregate species
T2 = groupsummary(T,{'sampling_event_id_res_check','day_of_year','lat','lon','all_obs_reported'},'sum',{'obs_count_0'});


Ts = groupsummary(T2,{'day_of_year','lat','lon'},'median',{'sum_obs_count_0'});
Ts=Ts(Ts.GroupCount>5,:);


%% Mapp it
lon = min(Ts.lon):dlatlon:max(Ts.lon);
lat = min(Ts.lat):dlatlon:max(Ts.lat);
M=nan(numel(lon),numel(lat),numel(unique(Ts.day_of_year)));
M(sub2ind( size(M),(Ts.lon-lon(1))/dlatlon+1,(Ts.lat-lat(1))/dlatlon+1,Ts.day_of_year-min(Ts.day_of_year)+1))=Ts.median_sum_obs_count_0;

figure;
tiledlayout('flow','TileSpacing','tight','Padding','tight')
for i_t=2:size(M,3)
    nexttile; 
    % tmp = M(:,:,i_t)'-nanmean(M,3)';
    tmp = log(M(:,:,i_t))';%-M(:,:,i_t-1)';
    imagesc(sort(unique(Ts.lon)),sort(unique(Ts.lat)),tmp,'AlphaData',~isnan(tmp))
    %Ts2=Ts(Ts.day_of_year==i_t,:);
    % scatter(Ts2.lon,Ts2.lat,500,Ts2.mean_sum_obs_count,'filled')
    % title(datestr(gext_dayY(i_t,i_y-1999)))
    axis tight equal; set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]); caxis([0 max(log(M(:)))]); %colorbar;
end

%%

%% Daily accumulation
Fs_accum=nan(g.nlat,g.nlon,366,21);
gext_dayY=NaT(366,21);
for i_y=2020:2020
    load(['data/flow/sim_' num2str(i_y) '.mat'],'Fd','gext')
    Fs_accum(:,:,:,i_y-1999) = mean(Fd.takingoff,4,'omitnan')+mean(Fd.landing,4,'omitnan');
    gext_dayY(:,i_y-1999) = gext.day;
end

%%
i_y = 2020;
i_tt = 256:280;

%% 
figure('position',[0 0 1000 600]); 
tiledlayout('flow','TileSpacing','tight','Padding','tight')
K = 1/(5^2)*ones(5);
for i_t=i_tt
    nexttile;
    tmp = conv2(-Fs_accum(:,:,i_t,i_y-1999)./g.area,K,'same');
    imagesc(g.lon,g.lat,tmp,'alphadata',~g.mask_water); 
    axis tight equal; set(gca,'ydir','normal'); borders('states','k'); axis([-125 -68 23 50]);colorbar; caxis([-300 300])
    title(datestr(gext_dayY(i_t,i_y-1999)))
end
