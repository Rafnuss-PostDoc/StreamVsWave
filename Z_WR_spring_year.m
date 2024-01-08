col2=brewermap([],'Paired');
clmap = brewermap([],'Spectral');
clmapd = crameri('berlin'); %hex2rgb(['#9e0142'; '#a00342'; '#a20543'; '#a40843'; '#a60a44'; '#a90d44'; '#ab0f45'; '#ad1245'; '#af1446'; '#b21746'; '#b41947'; '#b61c47'; '#b81e48'; '#bb2148'; '#bd2349'; '#bf2649'; '#c1284a'; '#c42b4b'; '#c62d4b'; '#c8304c'; '#ca324c'; '#cd354d'; '#cf374d'; '#d13a4e'; '#d33c4e'; '#d03c4d'; '#c73a4a'; '#bf3746'; '#b63543'; '#ad3240'; '#a5303d'; '#9c2d3a'; '#932b36'; '#8b2833'; '#822530'; '#79232d'; '#712029'; '#681e26'; '#5f1b23'; '#561920'; '#4e161d'; '#451419'; '#3c1116'; '#340f13'; '#2b0c10'; '#220a0c'; '#1a0709'; '#110506'; '#080203'; '#000000'; '#020507'; '#040b0f'; '#061017'; '#08161e'; '#0a1b26'; '#0c212e'; '#0e2635'; '#102c3d'; '#123145'; '#14374d'; '#163d54'; '#18425c'; '#1a4864'; '#1c4d6c'; '#1e5373'; '#20587b'; '#225e83'; '#24638a'; '#266992'; '#286f9a'; '#2a74a2'; '#2c7aa9'; '#2e7fb1'; '#3085b9'; '#3286bc'; '#3484bb'; '#3682ba'; '#387fb9'; '#3a7db8'; '#3b7bb6'; '#3d78b5'; '#3f76b4'; '#4174b3'; '#4371b2'; '#446fb1'; '#466db0'; '#486aaf'; '#4a68ae'; '#4c66ad'; '#4d63ab'; '#4f61aa'; '#515fa9'; '#535ca8'; '#555aa7'; '#5658a6'; '#5855a5'; '#5a53a4'; '#5c51a3'; '#5e4fa2']);

load('../BMM-US/data/density/inference-trans.mat')
load('../BMM-US/data/speed/inference-trans.mat')
addpath('../BMM-US/functions/')
vid = trans.f_inv(vidTS);

%% Load data

% define year
i_y=2021;

% Load flow data
load(['../BMM-US/data/flow/est_' num2str(i_y) '.mat'],'Fd','gext')

% load intra-night data
load(['../BMM-US/data/density/est_' num2str(i_y) '.mat'],'idt_y')
g_mw = repmat(~g.mask_water,1,1,numel(idt_y));
vy = nan(numel(g.lat), numel(g.lon), numel(idt_y)); vx=vy; rho = vy;

load(['../BMM-US/data/density/est_' num2str(i_y) '.mat'])
rho(g_mw) = EmbT; % bird/km^2 .* repmat(area(g.mask_water),1,g.nt); % bird
load(['../BMM-US/data/speed/est_uv_' num2str(i_y) '.mat'])
vx(g_mw) = Em.u/1000*60*60; % m/s -> km/h (+) east, (-) wes
vy(g_mw) = Em.v/1000*60*60; % m/s -> km/h (+) north, (-) south

gy.time = g.time(idt_y);
gy.day_id = g.day_id(idt_y)-min(g.day_id(idt_y))+1;
gy.day = g.day(unique(g.day_id(idt_y)));

vidt = vid(idt_y,:);
%% define period

t_tmp1= datetime('2021-4-24');%datetime('2015-10-7');
t_tmp2= datetime('2021-5-9');

id_day = find(gext.day<=t_tmp2 & gext.day>=t_tmp1);

%% Average daily for the period

id_day3 = find(gy.day<=t_tmp2 & gy.day>=t_tmp1);

rhod=nan(numel(g.lat), numel(g.lon), numel(id_day3));
vxd=rhod;
vyd=rhod;
for i = 1:numel(id_day3)
    id = gy.day_id==id_day3(i);
    rhod(:,:,i) = mean(rho(:,:,id),3,'omitnan');
    vxd(:,:,i) = mean(vx(:,:,id),3,'omitnan');
    vyd(:,:,i) = mean(vy(:,:,id),3,'omitnan');
end

%% Centroid 
tmp = padarray(Fd.landing,[1 1],0);
tmp(isnan(tmp))=0;
tmp = tmp + Fd.leaving;
tmp1 = padarray(Fd.takingoff,[1 1],0);
tmp1(isnan(tmp1)) = 0;
tmp1 = tmp1+Fd.entering;

weighted_centroid = nan(size(tmp,3),5);
for i=1:size(tmp,3)
    pl = regionprops(true(size(tmp(:,:,i))), abs(tmp(:,:,i)), 'WeightedCentroid');
    pt = regionprops(true(size(tmp1(:,:,i))), abs(tmp1(:,:,i)), 'WeightedCentroid');
    ptl = [pl.WeightedCentroid; pt.WeightedCentroid];
    ptl = [interp1(1:numel(g.lon), g.lon, ptl(:,1))  interp1(1:numel(g.lat), g.lat, ptl(:,2))];

    arclen = distance(ptl(1,2),ptl(1,1),ptl(2,2),ptl(2,1));
    weighted_centroid(i,:) = [deg2km(arclen) ptl(:)'];
end

%% Area
perc_area = .5;
tmp1 = sort(abs(reshape(Fd.landing,[],size(Fd.landing,3))));
tmp2 = sort(abs(reshape(Fd.takingoff,[],size(Fd.takingoff,3))));

id1 = cumsum(tmp1,'omitnan')>(perc_area*sum(tmp1,'omitnan'));
% TS.landing_area(id) =  sum(id1) .* g.dy * mean(g.dx);
tmp1(~id1)=nan;
thr1 = reshape(min(tmp1,[],1,'omitnan'),1,1,[]);

id2 = cumsum(tmp2,'omitnan')>(perc_area*sum(tmp2,'omitnan'));
% TS.takingoff_area(id) =  sum(id2) .* g.dy * mean(g.dx);
tmp2(~id2)=nan;
thr2 = reshape(min(tmp2,[],1,'omitnan'),1,1,[]);

landing_area = (abs(Fd.landing)-thr1)>0;
landing_area(repmat(g.mask_water,1,1,size(landing_area,3))) = 0;

takingoff_area = (abs(Fd.takingoff)-thr2)>0;
takingoff_area(repmat(g.mask_water,1,1,size(takingoff_area,3))) = 0;

% boundaries = bwboundaries(landing_area(:,:,100));


%% Daily

bd = load('borderdata.mat');
Anan = cellfun(@(x) [x(:);NaN],bd.lat(247:302),'un',0);
latb = cell2mat(Anan(:));
Anan = cellfun(@(x) [x(:);NaN],bd.lon(247:302),'un',0);
lonb = cell2mat(Anan(:));

figure('position',[0 0 1650 1200]); set(gcf, 'color', 'none');
tiledlayout(2,2,'TileSpacing','none','Padding','tight')

ax1=nexttile; hold on; xticks([]); yticks([]); box on;
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
tmp=nan(size(g.mask_water));
tmp(~g.mask_water) = 0;
im1=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
h2 = quivermc(g.LAT,g.LON,vxd(:,:,1),vyd(:,:,1),'linewidth',2,'density',10,'reference',100,'colormap',brewermap([],'YlOrBr'));
plot3m(latb,lonb,1000,"w");
t = textm(25,-120,datestr(gext.day(id_day(1))),'color','w','FontSize',40);

ax2=nexttile; hold on; xticks([]); yticks([]); box on;
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
im2=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
plot3m(latb,lonb,1000,"w"); colormap(ax2,clmapd)
textm(25,-120,"Change on ground",'color','w','FontSize',40);

ax3=nexttile; hold on; xticks([]); yticks([]); box on;
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
im3=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
plot3m(latb,lonb,1000,"w"); colormap(ax3,clmapd)
textm(25,-120,"Departure",'color','w','FontSize',40);


ax4=nexttile; hold on; xticks([]); yticks([]); box on;
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
im4=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
plot3m(latb,lonb,1000,"w"); colormap(ax4,clmapd)
textm(25,-120,"Landing",'color','w','FontSize',40);

clim(ax2,[-1 1]*3e5)
clim(ax3,[-1 1]*3e5)
clim(ax4,[-1 1]*3e5)

v = VideoWriter('figures/daily_2021.mp4','MPEG-4');
v.Quality = 50;
v.FrameRate = 2;
open(v)

for i_v=1:numel(id_day)
    delete(h2)
    axes(ax1)
    h2 = quivermc(g.LAT,g.LON,vxd(:,:,i_v),vyd(:,:,i_v),'linewidth',2,'density',10,'reference',60,'color',brewermap(1,'YlOrBr'));%'colormap',brewermap([],'YlOrBr'));

    im1.CData = smooth2a(rhod(:,:,i_v),1);
    im2.CData = smooth2a(Fd.takingoff(:,:,id_day(i_v)) + Fd.landing(:,:,id_day(i_v)),1);
    im3.CData = smooth2a(Fd.takingoff(:,:,id_day(i_v)),1);
    im4.CData = smooth2a(Fd.landing(:,:,id_day(i_v)),1);


    t.String=datestr(gext.day(id_day(i_v)));

    refreshdata
    drawnow
    % exportgraphics(gcf, "figures_paper/daily_maps_"+t.String+".png",'BackgroundColor','k')
    writeVideo(v,getframe(fig))
    % exportgraphics(gcf, "figures_paper/figure_1.png",'BackgroundColor','k')
end

% close(v)



%% Figure paper

i_v=3;


figure('position',[0 0 1650 1200]); hold on; xticks([]); yticks([]); box on;
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
setm(gca,'frame','off'); set(gca,'Visible','off') % set(gca, 'color', 'none');


% Export map
% plot3m(latb,lonb,1000,"k");
% set(gcf,'renderer','Painters')
% exportgraphics(gcf, "figures_paper/figure_1_map.eps")

clim([-1 1]*3e5);  colormap(clmapd)
%A=smooth2a(Fd.landing(:,:,id_day(i_v)),1);
% A=smooth2a(Fd.takingoff(:,:,id_day(i_v)),1);
%A=smooth2a(Fd.landing(:,:,id_day(i_v))+Fd.takingoff(:,:,id_day(i_v)),1);

A = smooth2a(rhod(:,:,i_v),1);
   
A=fillmissing(A,'nearest');
A=fillmissing(A,'nearest',2);
geoshow(g.LAT,g.LON,A,'DisplayType','surface');

h2 = quivermc(g.LAT,g.LON,vxd(:,:,i_v),vyd(:,:,i_v),'linewidth',2,'density',10,'reference',60,'color',brewermap(1,'YlOrBr'));%'colormap',brewermap([],'YlOrBr'));

% exportgraphics(gcf,"figures_paper/untitled5.eps","Resolution",300)




%%
fig=figure('position',[0 0 1650 1000]); % set(gcf, 'color', 'none');
tiledlayout(2,ceil(numel(id_day)/2),'TileSpacing','none','Padding','tight')
colormap(clmapd)
for i_v=1:(numel(id_day))
    nexttile; hold on; xticks([]); yticks([]); box on;
    % set(gca, 'color', 'none'); 
    set(gca,'Visible','off')
    tmp1=(-Fd.takingoff(:,:,id_day(i_v)));
    tmp2= (Fd.landing(:,:,id_day(i_v))); %./double(g.area)
    plot(-nansum(tmp1,2),g.lat,'color',clmapd(1,:),LineWidth=1.5)
    plot(nansum(tmp2,2),g.lat,'color',clmapd(end,:),LineWidth=1.5)
    d = nansum(tmp2,2)-nansum(tmp1,2);
    %plot(d,g.lat,'--k',LineWidth=2.5)
    d2=d;d(d<0)=0;
    plot(d2,g.lat,'-r',LineWidth=2.5)
    d2=d;d(d>0)=0;
    plot(d2,g.lat,'-b',LineWidth=2.5)
    text(0,g.lat(end)+1,datestr(gext.day(id_day(i_v)),'dd-mmm'),'color','k','FontSize',18,'HorizontalAlignment','center');
    xlim([-1 1]*2.5e7);
    xline(0)
    ylim([25 52])
end
% exportgraphics(gcf, "figures_paper/figure_2b.png")

 



















%%

fig=figure('position',[0 0 1650 1200]);  set(gcf, 'color', 'none');
tiledlayout('flow','TileSpacing','none','Padding','tight')
colormap(clmapd)
for i_v=1:(numel(id_day))
    nexttile; hold on; xticks([]); yticks([]); box on;
    axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
    mlabel off; plabel off; gridm off; framem; tightmap; box off;
    set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
    im2=geoshow(g.LAT,g.LON,smooth2a(Fd.takingoff(:,:,id_day(i_v)) + Fd.landing(:,:,id_day(i_v)),1),'DisplayType','surface');
    plot3m(latb,lonb,1000000000,"w",LineWidth=.1);

    tmp = weighted_centroid(id_day(i_v),:);
    plot3m(tmp(4:5),tmp(2:3),1000000001,'w-', linewidth=3)
    plot3m(tmp(4),tmp(2),10000000000,'o',MarkerFaceColor=clmapd(1,:), MarkerEdgeColor="w", Linewidth=2, MarkerSize=12)
    plot3m(tmp(5),tmp(3),10000000000,'o',MarkerFaceColor=clmapd(end,:), MarkerEdgeColor="w", Linewidth=2, MarkerSize=12)

    textm(25,-120,string(gext.day(id_day(i_v)),"dd-MMM"),'color','w','FontSize',24);
    clim([-1 1]*3e5)
end
exportgraphics(gcf, "figures/daily_maps_together.png",'BackgroundColor','k')



%% video

fig=figure('position',[0 0 1600 900]); tiledlayout('flow','TileSpacing','none','Padding','none');
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gcf, 'color', 'none'); set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
tmp=nan(size(g.mask_water));
tmp(~g.mask_water) = 0;
im=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
s=scatterm(radar.lat,radar.lon,100,vidt(id_day(1),:),'filled','MarkerEdgeColor','k');
clim([0 .4*max(rho(:,:, ismember(gy.day_id,id_day)),[],[1 2 3])])
bd = load('borderdata.mat');
Anan = cellfun(@(x) [x(:);NaN],bd.lat(247:302),'un',0);
latb = cell2mat(Anan(:));
Anan = cellfun(@(x) [x(:);NaN],bd.lon(247:302),'un',0);
lonb = cell2mat(Anan(:));
plot3m(latb,lonb,1000,"w");
h2 = quivermc(g.LAT,g.LON,vx(:,:,id_day(1)),vy(:,:,id_day(1)),'linewidth',2,'density',10,'reference',100,'colormap',brewermap([],'YlOrBr'));

t = textm(25,-120,datestr(gy.time(id_day(1))),'color','w','FontSize',40);

% p = plotm(nan,nan,"-ow",linewidth=2);
% p2 = plotm(nan,nan,"-or",linewidth=2);

v = VideoWriter('figures/spring_2021.mp4','MPEG-4');
v.Quality = 50;
v.FrameRate = 10;
open(v)

for i=2:numel(id_day)
    id = find(gy.day_id==id_day(i));
    cg_lat=nan(1,numel(id));
    cg_lon=nan(1,numel(id));
    cg_lat2=nan(1,numel(id));
    cg_lon2=nan(1,numel(id));
    for i_v = 1:numel(id)

        t.String=datestr(gy.time(id(i_v)),'dd-mmm-yyyy HH:MM');
        im.CData = smooth2a(rho(:,:,id(i_v)),1);

        delete(h2)
        tmp1 = vx(:,:,id(i_v));
        tmp1(isnan(rho(:,:,id(i_v))))=nan;
        tmp2 = vy(:,:,id(i_v));
        tmp2(isnan(rho(:,:,id(i_v))))=nan;
        try
            h2 = quivermc(g.LAT,g.LON,tmp1,tmp2,'linewidth',2,'density',10,'reference',60,'colormap',brewermap([],'YlOrBr'));
        end

        delete(s(:))
        id_nan=~isnan(vidt(id(i_v),:));
        s=scatterm(radar.lat(id_nan),radar.lon(id_nan),100,vidt(id(i_v),id_nan),'filled','MarkerEdgeColor','k');
% 
%         delete(p)
%         delete(p2)
%         if sum(isnan(rho(:,:,id(i_v))),'all')<=11644
%             tmp = rho(:,:,id(i_v)); tmp(isnan(tmp))=0;
%             props = regionprops(true(size(tmp)), tmp, 'WeightedCentroid');
%             cg_lon(i_v) = double(g.lon(round(props.WeightedCentroid(1))));
%             cg_lat(i_v) = double(g.lat(round(props.WeightedCentroid(2))));
%             p = plotm(cg_lat,cg_lon,"-k",linewidth=2);
% 
%             cg_lat2(i_v) = nansum(radar.lat' .* vidt(id(i_v),:)) ./ nansum(vidt(id(i_v),:));
%             cg_lon2(i_v) = nansum(radar.lon' .* vidt(id(i_v),:)) ./ nansum(vidt(id(i_v),:));
%             p2 = plotm(cg_lat2,cg_lon2,"-r",linewidth=2);
%         end
        refreshdata
        drawnow
        % pause(.1);

        writeVideo(v,getframe(fig))

        %[imind,cm] = frame2im(getframe(gcf));
        %   % Write to the GIF File
        %   if i_v == i_v1
        %       imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        %   else
        %       imwrite(imind,cm,filename,'gif','WriteMode','append');
        %   end
        %export_fig([filename num2str(i_v) '.png'], '-transparent')

    end
end

close(v)

%%
% Need to reun sinkSource. function for 2015 to run this script
simest="est";
saveit=false;
% edit sinksource


fig=figure('position',[0 0 1650 1200]);  set(gcf, 'color', 'none');
tiledlayout(11,2,'TileSpacing','none','Padding','tight')

nexttile([1 2]);t = text(0.5,0.5,"rest",'color','w','FontSize',40,'VerticalAlignment','middle','HorizontalAlignment','center');
set(gca,'color', 'none'); xticks([]); yticks([]); set(gca,'Visible','off')

ax1=nexttile([5 2]); hold on; xticks([]); yticks([]);
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
tmp=nan(size(g.mask_water));
tmp(~g.mask_water) = 0;
im1=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
h2 = quivermc(g.LAT,g.LON,vx(:,:,1),vy(:,:,1),'linewidth',2,'density',7,'reference',100,'colormap',brewermap([],'YlOrBr'));
plot3m(latb,lonb,1000,"w");
% t = textm(25,-120,datestr(gext.time(id_day(1))),'color','w','FontSize',40);

ax3=nexttile([5 1]); hold on; xticks([]); yticks([]); box on;
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
im3=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
plot3m(latb,lonb,1000,"w"); colormap(ax3,clmapd)
textm(25,-120,"Departure",'color','w','FontSize',40);


ax4=nexttile([5 1]); hold on; xticks([]); yticks([]); box on;
axesm('mercator','MapLatLimit',[23 50],'MapLonLimit',[-125 -67],'frame','off'); hold on;
mlabel off; plabel off; gridm off; framem; tightmap; box off;
set(gca, 'color', 'none'); setm(gca,'frame','off'); set(gca,'Visible','off')
im4=geoshow(g.LAT,g.LON,tmp,'DisplayType','surface');
plot3m(latb,lonb,1000,"w"); colormap(ax4,clmapd)
textm(25,-120,"Landing",'color','w','FontSize',40);

clim(ax1,[0 1]*400)
clim(ax3,[-1 1]*5e5)
clim(ax4,[-1 1]*5e5)

v = VideoWriter('figures/cum_2021.mp4','MPEG-4');
v.Quality = 50;
v.FrameRate = 10;
open(v)

vx(isnan(rho))=nan;
vy(isnan(rho))=nan;


for i=t_tmp1:t_tmp2
    [~,tmp] = min(abs(gext.time-i));
    id_day = find(gext.day_id==gext.day_id(tmp));


    for i_v=1:numel(id_day)
        delete(h2)
        axes(ax1)
        h2 = quivermc(g.LAT,g.LON,vx(:,:,id_day(i_v)),vy(:,:,id_day(i_v)),'linewidth',2,'density',7,'reference',60,'color',brewermap(1,'YlOrBr'));%'colormap',brewermap([],'YlOrBr'));

        im1.CData = smooth2a(rho(:,:,id_day(i_v)),1);
        im3.CData = smooth2a(-double(sum(takingoff(:,:,id_day(1:i_v)),3)),1);
        im4.CData = smooth2a(-double(sum(landing(:,:,id_day(1:i_v)),3)),1);

        t.String=datestr(gext.time(id_day(i_v)),'dd-mmm-yyyy HH:MM');

        refreshdata
        drawnow
        %pause(2)
        %exportgraphics(gcf, "figures/presentation/daily_maps_"+t.String+".png",'BackgroundColor','k')
        writeVideo(v,getframe(fig))
    end
end

close(v)






%%
Wlat = squeeze(sum(W,2,'omitnan'));

% full_night = squeeze(sum(rho>0,[1 2]))>12000;
sc = 2e6;%log(5e6);

figure('position',[0 0 1650 1200]); tiledlayout('flow','TileSpacing','none','Padding','tight')
for i_t = t_tmp1:t_tmp2
    ax1=nexttile;  hold on; box on;

    id_t = find(find(g.day==i_t)==gext.day_id);
    id_t=id_t(2:end);

    title(datestr(i_t,"dd-mmm"))
    % tmp = squeeze(nansum(takingoff(:,:,id_t),2));
    % tmp(tmp<10*mean(g.area,'all'))=0;
    imagesc(datenum(gext.time(id_t)),g.lat,Wlat(:,id_t))


    [~,id] = max(Wlat(:,id_t),[],'all');
    [lat_max, t_max] = ind2sub(size(Wlat(:,id_t)),id);
    plot(datenum(gext.time(id_t)), g.lat(lat_max)+((1:numel(id_t))-t_max) * 38 * 0.25 / 111,'-w')
    
    
   % tmp = squeeze(nansum(landing(:,:,id_t),2));
    %imagesc(datenum(gext.time(id_t)),g.lat,-tmp,'alphadata', -tmp./(sc/10))

    set(gca,"ydir","normal")
    clim([-1 1]*sc);
    % xlim(datenum([t_tmp1 t_tmp2]));
    datetick("x","HH:MM",'keeplimits'); axis tight
    ylim([27 47])
    %xticks([]); yticks([])

end
colormap(clmapd)
% exportgraphics(gcf, "figures_paper/nightly_landing_departure.png")













