col2=brewermap([],'Paired');
clmap = brewermap([],'Spectral');
clmapd = crameri('berlin'); %hex2rgb(['#9e0142'; '#a00342'; '#a20543'; '#a40843'; '#a60a44'; '#a90d44'; '#ab0f45'; '#ad1245'; '#af1446'; '#b21746'; '#b41947'; '#b61c47'; '#b81e48'; '#bb2148'; '#bd2349'; '#bf2649'; '#c1284a'; '#c42b4b'; '#c62d4b'; '#c8304c'; '#ca324c'; '#cd354d'; '#cf374d'; '#d13a4e'; '#d33c4e'; '#d03c4d'; '#c73a4a'; '#bf3746'; '#b63543'; '#ad3240'; '#a5303d'; '#9c2d3a'; '#932b36'; '#8b2833'; '#822530'; '#79232d'; '#712029'; '#681e26'; '#5f1b23'; '#561920'; '#4e161d'; '#451419'; '#3c1116'; '#340f13'; '#2b0c10'; '#220a0c'; '#1a0709'; '#110506'; '#080203'; '#000000'; '#020507'; '#040b0f'; '#061017'; '#08161e'; '#0a1b26'; '#0c212e'; '#0e2635'; '#102c3d'; '#123145'; '#14374d'; '#163d54'; '#18425c'; '#1a4864'; '#1c4d6c'; '#1e5373'; '#20587b'; '#225e83'; '#24638a'; '#266992'; '#286f9a'; '#2a74a2'; '#2c7aa9'; '#2e7fb1'; '#3085b9'; '#3286bc'; '#3484bb'; '#3682ba'; '#387fb9'; '#3a7db8'; '#3b7bb6'; '#3d78b5'; '#3f76b4'; '#4174b3'; '#4371b2'; '#446fb1'; '#466db0'; '#486aaf'; '#4a68ae'; '#4c66ad'; '#4d63ab'; '#4f61aa'; '#515fa9'; '#535ca8'; '#555aa7'; '#5658a6'; '#5855a5'; '#5a53a4'; '#5c51a3'; '#5e4fa2']);
addpath('../BMM-US/functions/')

load('../BMM-US/data/density/inference-trans.mat')

%% 
y = 2004:2021;

TS = table(g.day,VariableNames="day");
TS.year = year(TS.day);
TS.landing(:) = nan;
TS.takingoff(:) = nan;
TS.entering(:) = nan;
TS.leaving(:) = nan;

TS.NDTL(:) = nan;
LAT_NDTL = nan(numel(g.lat), 366, numel(y));

TS.landing_area(:) = nan;
TS.takingoff_area(:) = nan;
TS.overlap_area(:) = nan;
perc_area=[0:.01:.1 0.05:.05:1];
LD_AREA = nan(3, numel(perc_area), numel(y));

FD.landing=nan(g.nlat,g.nlon,numel(y),2);
FD.takingoff=FD.landing;

TS.wc_dist(:) = nan;

s=[2 2];
kernel=nan(numel(g.lat),numel(g.lon),height(radar));
for i_r = 1:height(radar)
    kernel(:,:,i_r) =  exp(- sqrt(double((abs(g.LON-radar.lon(i_r))/s(1)).^2 + (abs(g.LAT-radar.lat(i_r))/s(2)).^2)));
    % A = sqrt(double((abs(g.LON-radar.lon(i_r))).^2 + (abs(g.LAT-radar.lat(i_r))).^2));
end
FDR = nan(height(TS), height(radar),3);

corrp1 = nan(6, numel(y),2);

for i_y=1:numel(y)
    load(['../BMM-US/data/flow/est_' num2str(y(i_y)) '.mat'])

    % Remove outliar
    idt_rm = any(gext.day == [datetime(2006,3,24) datetime(2006,4,[19 20 25 29 30]) datetime(2006,4,[11 30]) datetime(2006,4,29) datetime(2012,3,19) datetime(2012,9,8) datetime(2012,10,2)],2);

    Fd.landing(:,:,idt_rm) = nan;
    Fd.takingoff(:,:,idt_rm) = nan;
    Fd.entering(:,:,idt_rm) = nan;
    Fd.leaving(:,:,idt_rm) = nan;
    MVT_day(:,:,idt_rm) = nan;
    Ts.landing_day(idt_rm) = nan;
    Ts.takingoff_day(idt_rm) = nan;
    Ts.entering_day(idt_rm) = nan;
    Ts.entering_day(idt_rm) = nan;


    % Smoothing
    for i=1:size(Fd.landing,3)
        Fd.landing(:,:,i) = smooth2a(Fd.landing(:,:,i),1);
        Fd.takingoff(:,:,i)= smooth2a(Fd.takingoff(:,:,i),1);
        MVT_day(:,:,i)= smooth2a(MVT_day(:,:,i),1);
    end

    % Find location in table
    [~,idt] = ismember(gext.day,TS.day);

    TS.landing(idt) = Ts.landing_day;
    TS.takingoff(idt) = Ts.takingoff_day;
    TS.entering(idt) = Ts.entering_day;
    TS.leaving(idt) = Ts.leaving_day;


    % NDTL
    diff = abs(Fd.takingoff)-abs(Fd.landing);
    norm = (abs(Fd.landing)+abs(Fd.takingoff)) / 2;
    % rd = abs(diff)./(norm*2);
    TS.NDTL(idt) = sum(abs(diff),[1 2],'omitnan') ./ sum(norm*2,[1 2],'omitnan');
    LAT_NDTL(:,1:numel(idt),i_y) = sum(abs(diff),2,'omitnan') ./ sum(norm*2,2,'omitnan');
    % figure; box on; grid on; plot(TS.day(idt),TS.turnover(idt)); ylabel("|Landing+Takingoff|/|Landing|+|Takingoff|")


    % Compute area of departure, arrival and overlap
    for i_p = 1:numel(perc_area)
        tmp1 = sort(abs(reshape(Fd.landing,[],size(Fd.landing,3))));
        tmp2 = sort(abs(reshape(Fd.takingoff,[],size(Fd.takingoff,3))));

        id1 = cumsum(tmp1,'omitnan')>(perc_area(i_p)*sum(tmp1,'omitnan'));
        % TS.landing_area(id) =  sum(id1) .* g.dy * mean(g.dx);
        tmp1(~id1)=nan;
        thr1 = reshape(min(tmp1,[],1,'omitnan'),1,1,[]);
    
        id2 = cumsum(tmp2,'omitnan')>(perc_area(i_p)*sum(tmp2,'omitnan'));
        % TS.takingoff_area(id) =  sum(id2) .* g.dy * mean(g.dx);
        tmp2(~id2)=nan;
        thr2 = reshape(min(tmp2,[],1,'omitnan'),1,1,[]);
        
        TS.landing_area(idt) =  sum(((abs(Fd.landing)-thr1)>0) .* g.area,[1 2]);
        TS.takingoff_area(idt) =  sum(((abs(Fd.takingoff)-thr2)>0) .* g.area,[1 2]);
        TS.overlap_area(idt) =  sum(((abs(Fd.takingoff)-thr2)>0 & (abs(Fd.landing)-thr1)>0) .* g.area,[1 2]);

        LD_AREA(1,i_p,i_y) = sum( TS.landing_area(idt) .* TS.landing(idt)) / sum(TS.landing(idt));
        LD_AREA(2,i_p,i_y) = sum( TS.takingoff_area(idt) .* TS.takingoff(idt)) / sum(TS.takingoff(idt));
        LD_AREA(3,i_p,i_y) = sum( TS.overlap_area(idt) .* TS.takingoff(idt)) / sum(TS.takingoff(idt));
    end

    % Seasonal map
    id_s = month(gext.day)>7;
    FD.landing(:,:,i_y,1) = nansum(Fd.landing(:,:,id_s),3);
    FD.takingoff(:,:,i_y,1) = nansum(Fd.takingoff(:,:,id_s),3);
    FD.landing(:,:,i_y,2) = nansum(Fd.landing(:,:,~id_s),3);
    FD.takingoff(:,:,i_y,2) = nansum(Fd.takingoff(:,:,~id_s),3);

    % Center of gravity
    tmp = padarray(Fd.landing,[1 1],0);
    tmp(isnan(tmp))=0;
    tmp = tmp + Fd.leaving;
    tmp1 = padarray(Fd.takingoff,[1 1],0);
    tmp1(isnan(tmp1)) = 0;
    tmp1 = tmp1+Fd.entering;

    for i=1:size(tmp,3)
        pl = regionprops(true(size(tmp(:,:,i))), abs(tmp(:,:,i)), 'WeightedCentroid');
        pt = regionprops(true(size(tmp1(:,:,i))), abs(tmp1(:,:,i)), 'WeightedCentroid');
        ptl = [pl.WeightedCentroid; pt.WeightedCentroid];
        ptl = [interp1(1:numel(g.lon), g.lon, ptl(:,1))  interp1(1:numel(g.lat), g.lat, ptl(:,2))];

        arclen = distance(ptl(1,2),ptl(1,1),ptl(2,2),ptl(2,1));
        TS.wc_dist(idt(i)) = deg2km(arclen);
    end

    % passage/change at radar location
     for i_r = 1:height(radar)
        kk = kernel(:,:,i_r);
        kk(isnan(Fd.takingoff(:,:,1))) = nan;
        kk = kk ./ sum(kk,"all","omitnan");
        
        FDR(idt,i_r,1) = squeeze(sum(Fd.takingoff./double(g.area) .* kk,[1 2],'omitnan')); % bird/km^2
        FDR(idt,i_r,2) = squeeze(sum(Fd.landing./double(g.area) .* kk,[1 2],'omitnan')); % bird/km^2
        FDR(idt,i_r,3) = squeeze(sum(MVT_day.*kk,[1 2],'omitnan')); % bird/km (summed over the night - averaged over spaced
     end

     % +1 day
     for i=0:5
        tmp1 = -reshape(Fd.landing(:,:,1:end-i),[],1);
        tmp2 = reshape(Fd.takingoff(:,:,(1+i):end),[],1);
        id=~isnan(tmp1)& ~isnan(tmp2);
        corrp1(i+1,i_y,1)=corr((tmp1(id)), (tmp2(id)));
        tmp2 = -reshape(Fd.landing(:,:,(1+i):end),[],1);
        id=~isnan(tmp1)& ~isnan(tmp2);
        corrp1(i+1,i_y,2)=corr((tmp1(id)), (tmp2(id)));
     end
     
end



%% Passage vs Chg on the ground

T = table(reshape(FDR(:,:,1),[],1), reshape(FDR(:,:,2),[],1), reshape(FDR(:,:,3),[],1),...
    reshape(repmat(TS.day,1,height(radar)),[],1),  reshape(repmat(string(radar.name)',numel(TS.day),1),[],1),...
    variablename=["takingoff", "landing", "passage", "day", "radar"]);
T.season = month(T.day)<7;
T(isnan(T.passage),:)=[];
T(T.passage<100,:)=[];
T.chg = -(T.landing+T.takingoff);
T.year = year(T.day);

x_x = logspace(0,log10(8000),100);
fitresults=fit(T.passage, T.chg,'poly1');
y_y_m=fitresults(x_x);
fitresults=fit(T.passage, T.chg.^2,'poly2');
y_y_std=sqrt(max(0,fitresults(x_x)));

a_axis = [0 8000 -1000 1000];
figure;tiledlayout(5,5,'TileSpacing','none','Padding','tight');
ax=nexttile([4,1]); 
histogram(T.chg,'EdgeColor','none','Orientation','horizontal');
ylim(a_axis(3:4)); ax.YGrid='on';% ax.XScale="log";
ylabel("Change on the ground (bird/km^2/night)")
nexttile([4,4]); hold on;
sz = max(T.passage/100,1);
id = datasample(1:numel(T.passage),20000,Replace=false, Weights=sz);
scatter(T.passage(id), T.chg(id), sz(id), 'o' ,'filled','markerfacealpha', 0.3); 
shadedErrorBar(x_x,y_y_m, y_y_std)
yline(0,'-k'); axis(a_axis);
box on; grid on; xticklabels(''); yticklabels('')
nexttile();
ax=nexttile([1,4]);
histogram(T.passage, 'EdgeColor','none')
xlim(a_axis(1:2));ax.XGrid='on';
xlabel("Passage (bird/km/night)")

% In terms of football pitch
% [1 2]./ (91.80*48.75/1000/1000) -> 223 447

% Secondary figures
figure; hold on;
plot(day(T.day,'doy'),T.takingoff,'.')
plot(day(T.day,'doy'),T.landing,'.')

quantile(T.takingoff,0.99)
median(groupsummary(T,"year","max","takingoff").max_takingoff)
median(groupsummary(T,"year","min","landing").min_landing)

figure; plot(reshape(FDR(:,:,1),[],1), reshape(FDR(:,:,2),[],1),'.k') 



%% Maps average around the radar location of take-off, landing, change and passage

Trsy = groupsummary(T,["radar", "season", "year"] ,"mean",["takingoff", "landing", "chg" "passage"]);
Trsy = join(Trsy,radar,LeftKeys="radar", RightKeys="name");

% Average over all years
% Mean year flux weighted by the number of day in the year available and
% then multiply by the max number of day to get total per season. 
Trs = groupsummary(Trsy,["radar", "season"] ,@(x,y) max(y) .* sum(x.*y)./sum(y) , {["mean_landing" "mean_takingoff" "mean_chg" "mean_passage"], ["GroupCount" "GroupCount" "GroupCount" "GroupCount"]});
Trs = join(Trs,radar,LeftKeys="radar", RightKeys="name");
Trs.takingoff = Trs.fun1_mean_takingoff_GroupCount;
Trs.landing = Trs.fun1_mean_landing_GroupCount;
Trs.chg = Trs.fun1_mean_chg_GroupCount;
Trs.passage = Trs.fun1_mean_passage_GroupCount;


figure('position',[0 0 900 900]); tiledlayout(4,2,'TileSpacing','tight','Padding','tight')
var = ["passage" "takingoff" "landing" "chg" ];
c_lim = [0 5e4; -1.5e4 1.5e4; -1.5e4 1.5e4; -2e3 2e3];
for i_v=1:numel(var)
    for i_s=1:-1:0
        f = fit([Trs.lon(Trs.season==i_s), Trs.lat(Trs.season==i_s)],Trs.(var(i_v))(Trs.season==i_s),'lowess');
        tmp = reshape(f([g.LON(:) g.LAT(:)]), size(g.LAT));
        tmp(g.mask_water)=nan;
        nexttile; hold on; 
        imagesc(g.lon, g.lat, tmp, alphadata=~isnan(tmp))
        axis equal tight; a_axis=axis; borders('states','k'); axis(a_axis);
        scatter(Trs.lon(Trs.season==i_s), Trs.lat(Trs.season==i_s), 100, Trs.(var(i_v))(Trs.season==i_s), 'filled','MarkerEdgeColor','k')
        clim(max(abs(tmp(:))).*[-1 1])
        title(var(i_v)+ " "+i_s); colorbar; box on;
        clim(c_lim(i_v,:))
        if i_v>1
            colormap(gca,clmapd)
        end
    end
end

%colormap(gca,clmapd)
%colormap(gca,clmap)

% Yearly change on the ground
var = ["mean_chg"];
for i_v=1:numel(var)
    for i_s=0:1
        figure('position',[0 0 900 600]); tiledlayout('flow','TileSpacing','tight','Padding','tight')
        for i_y=1:numel(y)
            Trs = Trsy(Trsy.year==y(i_y),:);
            f = fit([Trs.lon(Trs.season==i_s), Trs.lat(Trs.season==i_s)],Trs.(var(i_v))(Trs.season==i_s),'lowess');
            tmp = reshape(f([g.LON(:) g.LAT(:)]), size(g.LAT));
            tmp(g.mask_water)=nan;
            nexttile; hold on;
            imagesc(g.lon, g.lat, tmp, alphadata=~isnan(tmp))
            axis equal tight; a_axis=axis; borders('states','k'); axis(a_axis);
            scatter(Trs.lon(Trs.season==i_s), Trs.lat(Trs.season==i_s), 100, Trs.(var(i_v))(Trs.season==i_s), 'filled','MarkerEdgeColor','k')
            clim(max(abs(tmp(:))).*[-1 1])
            title(var(i_v)+ " "+i_s); colorbar; box on;
        end
    end
end

%% MAP
FD.landing(FD.landing==0)=nan;
FD.takingoff(FD.takingoff==0)=nan;

figure;  tiledlayout(2,2,'TileSpacing','tight','Padding','tight'); 
nexttile; imagesc(g.lon,g.lat, -smooth2a(mean(FD.landing(:,:,:,2)./double(g.area),3),1) ); set(gca,"YDir","normal"); title("Landing Spring")
nexttile; imagesc(g.lon,g.lat, -smooth2a(mean(FD.landing(:,:,:,1)./double(g.area),3),1) ); set(gca,"YDir","normal"); title("Landing Autumn")
nexttile; imagesc(g.lon,g.lat, smooth2a(mean(FD.takingoff(:,:,:,2)./double(g.area),3),1) ); set(gca,"YDir","normal"); title("Takingoff Spring")
nexttile; imagesc(g.lon,g.lat, smooth2a(mean(FD.takingoff(:,:,:,1)./double(g.area),3),1) ); set(gca,"YDir","normal"); title("Takingoff Autumn")


%% Max departure
k=5;
Lmax=nan(numel(y),k,4);
for i_y=1:numel(y)
    Lmax(i_y,:,3)=maxk(abs(TS.landing(TS.year==y(i_y) & month(TS.day)<7)),k);
    Lmax(i_y,:,2)=maxk(abs(TS.landing(TS.year==y(i_y) & month(TS.day)>=7)),k);
    Lmax(i_y,:,4)=maxk(abs(TS.takingoff(TS.year==y(i_y) & month(TS.day)<7)),k);
    Lmax(i_y,:,1)=maxk(abs(TS.takingoff(TS.year==y(i_y) & month(TS.day)>=7)),k);
end
figure; hold on
boxplot(Lmax(:,:,1),'PlotStyle','compact')
boxplot(Lmax(:,:,2),'PlotStyle','compact')
boxplot(Lmax(:,:,3),'PlotStyle','compact')
boxplot(Lmax(:,:,4),'PlotStyle','compact')
ylim([0 2e9])

figure; hold on; box on; grid on;
bar(squeeze(quantile(Lmax,.5)))  
errorbar(1:5,squeeze(quantile(Lmax,.5)), squeeze(quantile(Lmax,.5))-squeeze(quantile(Lmax,.25)), squeeze(quantile(Lmax,.75))-squeeze(quantile(Lmax,.5))); 
xlabel("Top Nights for each year"); ylabel("Number of birds")
legend(["Landing Spring", "Landing Autumn", "Takingoff Spring", "Takingoff Autumn"])


%% Absolute Relative Difference
w = (-TS.landing+TS.takingoff)/2;

nansum(TS.NDTL.*w)./ nansum(w)

nanstd(TS.NDTL, w)

[histw, intervals] = histwc(TS.NDTL, w, 30);
figure; hold on;
histogram(TS.NDTL)
bar(intervals, histw,1)
xlim([0 1]); grid on; box on; yticks([])
xlabel("Normalized Absolute Difference")

%
h = fspecial3('gaussian',[10 10 1],4);
LAT_NDTL_sm = imfilter(LAT_NDTL,h);
tmp = mean(LAT_NDTL_sm,3,'omitnan');
tmp2 = std(LAT_NDTL_sm,[],3,'omitnan');
%tmp = smooth(mean(LAT_NDTL,3),2);

figure; hold on; box on; grid on;
id_lat = find(any(g.lat' == 33:5:45,2));
t = datenum(TS.day(TS.year==2000));
for i_lat=1:numel(id_lat)
    shadedErrorBar(t,tmp(id_lat(i_lat),:)', tmp2(id_lat(i_lat),:)')
end
legend(num2str(g.lat(id_lat)'))
datetick('x'); axis tight;
ax=gca; ax.XTick = datenum(datetime(2000,1:12,1));
ylabel("Aboslute Relative Difference")

% 
figure; plot((-TS.landing+TS.takingoff)/2, TS.NDTL,'.k')
set(gca,"xscale","log");
xlabel("(Departure + Arrival ) / 2")
ylabel("Normalized Absolute Difference")

%
figure; box on; grid on;
plot(day(TS.day,'doy'),(-TS.landing+TS.takingoff)/2,'.k')
ylabel("(Departure + Arrival ) / 2")

%% Area per night
figure; 
tmp = mean(LD_AREA,3,'omitnan')';
tmp = [tmp(:,1)-tmp(:,3) tmp(:,3) tmp(:,2)-tmp(:,3)];
area(1-perc_area, tmp)
newcolors = [clmapd(1,:); 0.85 0.85 0.85;clmapd(end,:)];
colororder(newcolors)
box on; grid on; 
% set(gca,"yscale","log")

% figure; plot(1-perc_area,tmp./sum(tmp,2))




%% Distance centroid
figure; box on; grid on;
id = ~isnan(TS.wc_dist);
scatter(day(TS.day(id),'doy'),TS.wc_dist(id),(-TS.landing(id)+TS.takingoff(id))/10^7,'ok','filled','MarkerFaceAlpha',.4)
ylabel("Centroid displacement (km)")

w = (-TS.landing+TS.takingoff)/2;
[histw, intervals] = histwc(TS.wc_dist(TS.wc_dist<=300), w(TS.wc_dist<=300),30);
figure('position',[0 0 600 200]); hold on; xticks([]);
barh(intervals, histw, 1)
yline(nansum(TS.wc_dist.*w./nansum(w)))
yline(nansum(TS.wc_dist(month(TS.day)<7).*w(month(TS.day)<7)./nansum(w(month(TS.day)<7))))
yline(nansum(TS.wc_dist(month(TS.day)>7).*w(month(TS.day)>7)./nansum(w(month(TS.day)>7))))
ylim([0 300]); grid on; box on; set(gca,"ydir","reverse")
ylabel("Distance centroid landing - take-off")

figure; plot(TS.wc_dist,TS.NDTL,'.')
xlabel("Centroid displacement (km)")
ylabel("NDTL")

mean(TS.wc_dist(w>quantile(w,.95)))

i=find(TS.day==datetime(2021,5,9));


%% Night duration
% 8h30-13h per night. 
tmp = datetime(2021,1,1):(1/24/60):datetime(2022,1,1);
NNT = twilightNNT(tmp, radar.lon, radar.lat);
NNT(NNT<thr_nnt(1)|NNT>thr_nnt(2))=nan;
tmp2 = table(NNT(:), reshape(repmat(tmp',1,height(radar)),[],1),  reshape(repmat(string(radar.name'),numel(tmp),1),[],1),...
    VariableNames=["nnt","time","radar"]);
tmp2.doy=day(tmp2.time,'doy');

tmp3 = groupsummary(tmp2,["doy", "radar"],@(x) mean(~isnan(x)),"nnt") ;

plot(tmp3.doy,tmp3.fun1_nnt*24,'.k'); grid on; box on; axis tight;
ylabel("Night duration"); xlabel("Day of year")


T.doy=day(T.day,"doy");
tmp4 = join(T,tmp3,LeftKeys=["doy", "radar"], RightKeys=["doy", "radar"]);


sum(tmp4.fun1_nnt .* tmp4.passage) ./ sum(tmp4.passage) * 24

id = ~isnan(vin) & ~isnan(vid) & vin>0;
% sum(vin(id) .* vid(id)) ./ sum(vid(id))*60*60/1000

%% +1 day

figure('position',[0 0 800 350]); hold on; grid on;
y = [median(corrp1(:,:,2),2) median(corrp1(:,:,1),2)];
ymin = y-[quantile(corrp1(:,:,2),.1,2) quantile(corrp1(:,:,1),.1,2)];
ymax = [quantile(corrp1(:,:,2),.9,2) quantile(corrp1(:,:,1),.9,2)]-y;
bar(y)
errorbar(1:6,y(:,1),ymin(:,1),ymax(:,1))
errorbar(1:6,y(:,2),ymin(:,2),ymax(:,2))

xticklabels("+"+string(0:5))
ylabel("Correlation"); xlabel("Landing at day 0 and Landing at day + x")