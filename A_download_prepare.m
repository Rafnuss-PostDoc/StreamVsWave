%%
% This script load the eBird data from the sql database and perform basic
% preliminary operation. It also load the equivalent WR dataset 

%% Loading data
path_bmmus="/Users/raphael/Library/CloudStorage/Box-Box/BMM-US";
load(path_bmmus+'/data/density/inference-trans.mat','g','radar')
addpath(genpath(path_bmmus+'/functions/'))
path_matlabebird="/Users/raphael/Library/CloudStorage/Box-Box/MatlabeBird";
addpath(genpath(path_matlabebird+'/functions/'))

% Reading nocturnal species

opts = delimitedTextImportOptions("NumVariables", 9);
opts.VariableNames = ["SPECIES_CODE", "SCI_NAME", "PRIMARY_COM_NAME", "group", "family", "population_size_2017", "BodyMassValue", "biomass", "DetectRadarNight"];
opts.VariableTypes = ["string", "string", "string", "categorical", "categorical", "double", "double", "double", "categorical"];

speNoc = readtable('data/SpeciesList-Wee_Hao_530.csv',opts);
speNoc=speNoc(speNoc.DetectRadarNight=="TRUE",:);
speNoc=speNoc(speNoc.group=="landbird",:);

% Definition of the cities

cities = readtable("data/cities.csv",'TextType','string');

% All city tested
% city_keep = ["Minneapolis", "Madison", "Ithaca", "Boston", "Chicago", "New York", "Philadelphia", "Cincinnati","Washington","Kansas City","St. Louis","Nashville-Davidson","Pittsburgh","Atlanta","Dallas","Detroit","Asheville","San Antonio","Gainesville","Denver","Houston","Austin", "Indianapolis","Ithaca (close)"];

% cities with enough data
city_keep = ["Minneapolis", "Madison", "Ithaca", "Boston","Detroit", "Chicago", "New York","Pittsburgh","Philadelphia", "Denver","Cincinnati","Kansas City","Washington","Asheville","Atlanta","Dallas","Austin","Houston"];

% filter city
cities = cities(ismember(cities.name,city_keep),:);

% Add ithaca missing
cities = [cities; {"Ithaca", 42.445478, -76.495112}];

% Sort by latitude
cities = sortrows(cities,"latitude","descend");









%% Query parameter definition
q_cond.longitude = [-1 100];
q_cond.longitude = [-1 100];
q_cond.effort_distance_km = [0 10];
q_cond.effort_hrs = [10/60 5];
q_cond.num_observers = [0 10000];
q_cond.day_of_year = [90 155];
q_cond.year = [2010 2050];
q_cond.cci = [-2 100];

radius = 0.75;

q_col_check=["checklist_id","observer_id","latitude","longitude","obs_dt","effort_hrs","effort_distance_km","num_observers","number_of_species","cci", "cds_tp", "cds_lcc", "cds_i10fg", "cds_t2m","cds_msl", "protocol_id"];
q_col_obs=["checklist_id","species_code","obs_count"];


%% Load data 
% Around 2.5 hr
for i_c = 1:height(cities)
    try
        q_cond.longitude = cities.longitude(i_c) + [-radius radius];
        q_cond.latitude = cities.latitude(i_c) + [-radius radius];
        [res_check, res_obs] = query_erd(q_cond,q_col_check=q_col_check);
      
        save("data/eBird/"+cities.name(i_c)+".mat",'res_check','res_obs','q_cond')
    catch
        disp("error with " + cities.name(i_c))
    end
end
 
%% Read and process.75
Tr=cell(height(cities),1);
T2r=cell(height(cities),1);

for i_c = 1:height(cities)
    disp("Loading: "+cities.name(i_c))
    load("data/eBird/"+cities.name(i_c)+".mat")

    % Filter distance in a radius
    res_check = res_check(sqrt((res_check.latitude-mean(q_cond.latitude)).^2 + (res_check.longitude-mean(q_cond.longitude)).^2 )<radius,:);
   
    T_obs_dt_utc = res_check.obs_dt; T_obs_dt_utc.TimeZone="UTC";
    [~,twlSet,twlRise,twDate] = twilightNNT(T_obs_dt_utc, mean(q_cond.longitude), mean(q_cond.latitude));
    twlSet.TimeZone = res_check.obs_dt.TimeZone;
    twlRise.TimeZone = res_check.obs_dt.TimeZone;
    twDate.TimeZone = res_check.obs_dt.TimeZone;
    Rise = interp1(dateshift(twlRise,"start",'day'),twlRise,dateshift(res_check.obs_dt,"start",'day'), 'nearest');
    Set = interp1(dateshift(twlSet,"start",'day'),twlSet,dateshift(res_check.obs_dt,"start",'day'), 'nearest');
    res_check.ss = min(hours(res_check.obs_dt-Rise),0) + max(hours(res_check.obs_dt-Set),0);
    res_check.hours_since_sunrise = minutes(res_check.obs_dt-Rise)/60;
    res_check.obs_dt.TimeZone = '';
    % figure; plot(T.obs_dt,T.ss,'.k')

    % Combine checklist and obs
    T = outerjoin(res_check,res_obs, Keys="checklist_id", MergeKeys=true, Type="left");
    T = sortrows(T,'obs_dt');
    T.name(:) = string(cities.name(i_c));
    T.num_observers=double(T.num_observers);
    T.obs_dt_day = dateshift(T.obs_dt,'start','day'); % Add day of time
    % height(T) 


    % Limit max count
    disp(round(mean(T.obs_count>100)*100,3)+"% of the obs data with a count above 100. count -> 100")
    T.obs_count(T.obs_count>100) = 100;

    
    % Check species list
%     Tsplist = groupsummary(T,{'species_code'},'sum',{'obs_count'});
%     [~,id]=ismember(Tsplist.species_code,speNoc.SPECIES_CODE);
%     Tsplist.RCS(id~=0) = speNoc.BodyMassValue(id(id~=0)).^(2/3);
%     Tsplist.w = round(Tsplist.RCS .* Tsplist.sum_obs_count ./ sum(Tsplist.RCS .* Tsplist.sum_obs_count)*100,2);
%     sortrows(Tsplist,'w','descend');
    
    % Filter by species by setting the obs_count to nan for all non
    % landbird nocturnal migrant
    % remove non-migratory species
    [~,id]=ismember(T.species_code,speNoc.SPECIES_CODE);
    disp(round(mean(id==0)*100)+"% of the obs data from an non-landbird nocturnal migrant. count -> nan")
    T.obs_count(id==0) = nan;
    % T.obs_count(T.species_code=="normoc")=nan;
    % T.obs_count(T.species_code=="amerob")=nan;
    % T.obs_count(T.species_code=="cedwax")=nan;
    % T.obs_count(T.species_code=="daejun")=nan;
    % T.obs_count(T.species_code=="whtspa")=nan;
    % T.obs_count(T.species_code=="sonspa")=nan;
    % remove non-warbler
    % warbler = {'yerwar','yelwar','ovenbi1','amered','chswar','foxspa','magwar','palwar','btnwar','tenwar','norwat','bkbwar','btbwar','bkbwar','leafly','buwwar','pinwar','bawwar'};
    % T.obs_count(~ismember(T.species_code,warbler))=nan;


    % Compute RCS for each species 
    T.obs_rcs(:)=nan;
    T.obs_rcs(id>0) = T.obs_count(id>0) .* speNoc.BodyMassValue(id(id>0)).^(2/3);
    

    % Account for zero: T.obs_count -> nan = not recorded, 0: X (presnce only)
    % set "X" to mean of the group
    % G=findgroups(T(:, ismember(T.Properties.VariableNames, {'dt','species_code'})));
    % obs_count_GROUP_Not0 = splitapply(@(x) quantile(x(x>0),q_cond.radius), T.obs_count, G);
    % T.obs_count_0 = T.obs_count;
    % T.obs_count_0(T.obs_count==0) = obs_count_GROUP_Not0(G(T.obs_count==0));
    list_event_id_presence_only = unique(T.checklist_id(T.obs_count==0));
    id = ismember(T.checklist_id,list_event_id_presence_only);
    disp(round(mean(id)*100)+"% of the data from a checklist not reporting all counts. DELETED")
    T(id,:)=[];


    % PART DIVERSITY
    Tr{i_c}=T;
    
    % PART SUM COUNT
    % Group observation in checklist
    T2 = groupsummary(T,[res_check.Properties.VariableNames(:)',{'hours_since_sunrise'},{'name'}],'sum',{'obs_rcs','obs_count'});

    % Assing to data structure
    T2r{i_c}=T2;

    % PART SPECIES SUM
    % T3 = groupsummary(T,["obs_dt_day", "species_code"], "sum", "obs_count");
end

save("data/eBird/T2r.mat",'T2r')
save("data/eBird/Tr.mat",'Tr','-v7.3')






























%% Radar data
load(path_bmmus+'/data/density/inference-trans.mat','g','radar')
s=[2 2];
kernel=nan(size(g.LON,1),size(g.LON,2),height(cities));
for i_c = 1:height(cities)
    A = exp(- sqrt(double((abs(g.LON-mean(cities.longitude(i_c)))/s(1)).^2 + (abs(g.LAT-mean(cities.latitude(i_c)))/s(2)).^2)));
    % B = sqrt((g.LON-mean(cities.longitude(i_c))).^2 + (g.LAT-mean(cities.latitude(i_c))).^2)<=radius;
    % A(g.mask_water)=0;
    %B(g.mask_water)=0;
    %A(B) = sum(A(~B))/sum(B,'all');
    % sum(A(B))./sum(A,'all')
    kernel(:,:,i_c) = A;
end


ts = datetime(2010,1,1):datetime(2022,1,1);
Fd_takingoff = nan(numel(ts),height(cities));
Fd_landing = nan(numel(ts),height(cities));
Fd_diff = nan(numel(ts),height(cities));
Fd_mvt = nan(numel(ts),height(cities));
for i_y=2010:2021
    i_y
    load(['../BMM-US/data/flow/est_' num2str(i_y) '.mat'])
    
    for i_c = 1:height(cities)
        kk = kernel(:,:,i_c);
        kk(isnan(Fd.takingoff(:,:,1))) = nan;
        kk = kk ./ sum(kk,"all","omitnan");
        
        idt = find(ts==datetime(i_y,1,1))+(0:size(Fd.takingoff,3)-1);
        
        Fd_takingoff(idt,i_c)=squeeze(sum(Fd.takingoff./double(g.area) .* kk,[1 2],'omitnan')); % bird/km^2
        Fd_landing(idt,i_c)=squeeze(sum(Fd.landing./double(g.area) .* kk,[1 2],'omitnan')); % bird/km^2
        Fd_diff(idt,i_c)=cumsum(-Fd_takingoff(idt,i_c)-Fd_landing(idt,i_c));

        Fd_mvt(idt,i_c) =  squeeze(sum(MVT_day.*kk,[1 2],'omitnan')); % bird/km (summed over the night - averaged over spaced
    end
end

% save
save("data/radar_cities","cities","ts","Fd_diff","Fd_takingoff","Fd_landing","Fd_mvt")





























return

%% Create figure for paper.

% eBird coverage map
path_bmmus="/Users/rafnuss/Library/CloudStorage/Box-Box/BMM-US";
load(path_bmmus+'/data/density/inference-trans.mat','g','radar')

% Load all 
if false
    q_cond.longitude = [-125 -66];
    q_cond.latitude = [25 50];
    q_cond.day_of_year = [90 150];
    q_col_check=["checklist_id","observer_id","latitude","longitude","obs_dt","effort_hrs","effort_distance_km","num_observers","number_of_species","cci"];
    res_check = query_erd(q_cond,q_col_check=q_col_check);
    % save('data/eBird/entireUS.mat','res_check','q_cond','q_col_check');
else
    load('data/eBird/entireUS.mat')
end

TeUS=res_check;

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

M(round(A.id)) = A.GroupCount/(dlatlon*111)^2; 

figure('position',[0 0 900 500]);  set(gcf, 'color', 'none');
tiledlayout(1,1,'TileSpacing','none','Padding','tight')
hold on; xticks([]); yticks([]);
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off; 
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
M(M==0)=nan;
im1=geoshow(LAT,LON,(M'),'DisplayType','texturemap');
colormap(copper)

bordersm('states','color',.3.*[1 1 1]); 
cities_all = readtable("data/cities.csv",'TextType','string');
city_keep = ["Minneapolis", "Madison", "Ithaca", "Boston", "Chicago", "New York", "Philadelphia", "Cincinnati","Washington","Kansas City","St. Louis","Nashville-Davidson","Pittsburgh","Atlanta","Dallas","Detroit","Asheville","San Antonio","Gainesville","Denver","Houston","Austin", "Indianapolis","Ithaca (close)"];
cities_all = cities_all(ismember(cities_all.name,city_keep),:);
circlem(cities_all.latitude,cities_all.longitude,111*0.75,'edgecolor',[255 151 0]./255,'linewidth',1)
circlem(cities.latitude,cities.longitude,111*0.75,'edgecolor',[0.8080    0.1775    0.0166],'linewidth',2)
textm(cities.latitude+1.2,cities.longitude,cities.name,'HorizontalAlignment','center','color',[0.8080    0.1775    0.0166],'fontsize',12)

scatterm(radar.lat, radar.lon,80,'.y')
caxis([0 20])
c=colorbar; c.Color="w";

% exportgraphics(gcf, "figures_paper/ebird_map.png",'BackgroundColor','k')


%% 
figure('position',[0 0 1650 1200]);  set(gcf, 'color', 'none');
tiledlayout(1,1,'TileSpacing','none','Padding','tight')
hold on; xticks([]); yticks([]);
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off; 
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')

im1=geoshow(g.LAT,g.LON,kernel(:,:,3),'DisplayType','texturemap');
circlem(42.4553, -76.4772,111*q_cond.radius,'edgecolor',[0.8080    0.1775    0.0166],'linewidth',2)

bordersm('states','w'); 

% exportgraphics(gcf, "figures_paper/WR_map.png",'BackgroundColor','k')

%% Spatial change of distribution

dlatlon=.02;
figure(Position=[0,0,1400,800]); tiledlayout('flow','TileSpacing','tight','Padding','tight')
for i_c=1:height(cities)
    nexttile;
    
    T2c = T2(T2.name==cities.name(i_c) ,:); % & year(T2.obs_dt)>=2015 weekend
    T2c.lat = round(T2c.latitude/dlatlon)*dlatlon;
    T2c.lon = round(T2c.longitude/dlatlon)*dlatlon;
    
    A=groupcounts(T2c,{'lat','lon'});
    lon=min(A.lon):dlatlon:max(A.lon);
    lat=min(A.lat):dlatlon:max(A.lat);
    [LAT,LON]=meshgrid(lat,lon);
    M = nan(numel(lon), numel(lat));
    A.lon_id = (A.lon-lon(1)) / dlatlon+1;
    A.lat_id = (A.lat-lat(1)) / dlatlon+1;
    A.id = sub2ind(size(M),round(A.lon_id),round(A.lat_id));
    
    M(round(A.id)) = A.GroupCount; 
    M = M ./ sum(M,'all','omitnan');
    imagesc(lon,lat,M')
    title(cities.name(i_c))
    axis equal tight square; a_xis=axis;borders('states','w'); axis(a_xis)
    set(gca,ydir="Normal", XTick=[], YTick=[])
    clim([0 .05]);
end



%% Treemap 
% load("data/eBird/Tr.mat")
T = vertcat(Tr{:});

[~,id]=ismember(T.species_code,speNoc.SPECIES_CODE);
disp(round(mean(id==0)*100)+"% of the obs data from an non-landbird nocturnal migrant. count -> nan")
   
Tsplist = groupsummary(T,{'species_code'},'sum',{'obs_count'});
Tsplist = sortrows(Tsplist(Tsplist.sum_obs_count>0,:),'sum_obs_count','descend');

Tsplist.label = Tsplist.species_code;
tmp = Tsplist.sum_obs_count./sum(Tsplist.sum_obs_count);
tmp2 = cumsum(tmp)
Tsplist.label(tmp2>.86)="";
Tsplist.label(tmp>0.01)=Tsplist.label(tmp>0.01)+newline+" "+ round(tmp(tmp>0.01)*100)+"%";
figure('position',[0 0 800 450]);
rec=treemap(Tsplist.sum_obs_count,2,1);
plotRectangles(rec, Tsplist.label)
outline(rec)