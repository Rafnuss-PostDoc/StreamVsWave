
% run B_    
run("B_construct_daily_table.m")

%% Same obs-location

T2 = sortrows(T2,{'longitude','latitude','observer_id','obs_dt_day'});

% Not dealing with same-day checklist well.
id = (~diff(T2.longitude) & ~diff(T2.latitude) & ~diff(T2.observer_id) & days(diff(T2.obs_dt_day))==1 );

if false
% filtering
    val = ["effort_hrs", "effort_distance_km", "hours_since_sunrise"];
    thr = [30/60 .4 3];
    dc = abs(diff(T2.cp));
    
    figure; tiledlayout('flow')
    for i=1:numel(val)
        nexttile;
        tmp = abs(diff(T2.(val(i))));
        histogram(tmp(id),'Normalization','cdf')
        xline(thr(i),"r",LineWidth=3)
        xlabel("delta "+val(i))
        nexttile;
        plot(tmp(id),dc(id),'.k');
        xlabel("delta "+val(i)); ylabel("delta number of bird")
        l=lsline; l.Color="r"; l.LineWidth=2;
        %ylim([0 50])
    end
    
    id2 = id & abs(diff(T2.effort_hrs))<30/60;
    sum(id2)/sum(id)
    id2 = id & abs(diff(T2.effort_distance_km))<.4;
    sum(id2)/sum(id)
    id2 = id & abs(diff(T2.hours_since_sunrise)) < 3;
    sum(id2)/sum(id)
end


T20 = T2(find(id),:);
T21 = T2(find(id)+1,:);

T2_diff = table(T20.obs_dt_day, T20.cci, T21.cp-T20.cp, T21.effort_distance_km-T20.effort_distance_km, T21.effort_hrs-T20.effort_hrs,T21.hours_since_sunrise-T20.hours_since_sunrise, T20.checklist_id, T21.checklist_id,T21.name,...
    'VariableNames',{'ts','cci','dcp','ddist','dhrs','dhours_since_sunrise','checklist_id_0','checklist_id_1','name'} );







%% For number only

% Combine with radar
[G,TID]=findgroups(T2_diff(:,["ts","name"]));
TID.ebird_m=splitapply(@median,T2_diff.dcp,G);
TID.ebird_n=splitapply(@numel,T2_diff.dcp,G);
TID.ebird_q25=splitapply(@(x) quantile(x,.25),T2_diff.dcp,G);
TID.ebird_q75=splitapply(@(x) quantile(x,.75),T2_diff.dcp,G);
TID.ebird_d=splitapply(@(x) {x} ,T2_diff.dcp,G);

Tk_diff = outerjoin(TID,Trd,MergeKeys=true,RightKeys=["name","obs_dt_day"],LeftKeys=["name","ts"]);
Tk_diff = renamevars(Tk_diff,"ts_obs_dt_day","ts");




%% Figure
for i_c=3;%height(cities)
figure('position',[0 0 1650 850]); tiledlayout('flow','TileSpacing','tight','Padding','tight'); set(gcf, 'color', 'k');
 y_eb_max = 10*nanmean(abs((Tk_diff.ebird_m(Tk_diff.name==cities.name(i_c)))));
    y_wr_max = y_eb_max*8;
    for i_y=2016:2021
        Tkt = Tk_diff(year(Tk_diff.ts)==i_y&Tk_diff.name==cities.name(i_c)&Tk_diff.ebird_n>0,:);
           
        nexttile;hold on; title(i_y,'color','w')
        set(gca, 'color', 'k');
        set(gca,'XColor','w');
        set(gca,'yColor','w');
        
        %
        %bar(ts, Fd_dens(:,i_c), 1,'FaceColor',[.5 .5 .5]);
        bar(Tkt.ts,-Tkt.takingoff,1,'FaceColor',clmap(end,:),'FaceAlpha',.5)
        bar(Tkt.ts,-Tkt.landing,1,'FaceColor',clmap(1,:),'FaceAlpha',.5)
        bar(Tkt.ts,-Tkt.takingoff-Tkt.landing,'w')
           
        ylim([-y_wr_max y_wr_max])
        % xticklabels(''); yticklabels('');
        datetick('x','mmm','keepticks')
        box on; grid on;
        yyaxis right; 
    
        a = cumsum(Tkt.ebird_n,'omitnan');
        b = zeros(1,a(end));
        b(a - Tkt.ebird_n +1)=1;
        scatter(Tkt.ts(cumsum(b)), vertcat(Tkt.ebird_d{:}),1,'y','filled')
        errorbar(Tkt.ts, Tkt.ebird_m, Tkt.ebird_m-Tkt.ebird_q25, Tkt.ebird_q75-Tkt.ebird_m,'y',"LineStyle","none")
        scatter(Tkt.ts, Tkt.ebird_m, Tkt.ebird_n*4,'y','filled')
        yline(0,'--y',LineWidth=2)
        ylim([-y_eb_max y_eb_max])
        set(gca,'yColor','w');
        % yticklabels(''); 
    end
 % exportgraphics(gcf, "figures/same_obs_cumulative_"+num2str(i_c)+"_"+cities.name(i_c)+"_autumn.png",'BackgroundColor','k')
end


%% Daily correlation
thr_n=50;
c = nan(height(cities),3);

for i_c=1:height(cities)
    Tkt=Tk_diff;
    Tkt = Tkt(Tkt.ebird_n>thr_n&Tkt.name==cities.name(i_c),:);

    if sum(Tkt.ebird_n)>0
        c(i_c,3) = sum(Tkt.ebird_n);
        c(i_c,1) = corrW(Tkt.mvt,abs(Tkt.ebird_m),Tkt.ebird_n);
        c(i_c,2) = corrW(Tkt.takingoff+Tkt.landing, Tkt.ebird_m, Tkt.ebird_n);
    end
end

% 
if false
    y_lim=8;
    for i_c=3%1:height(cities)
        figure('position',[0 0 850 850]);  tiledlayout(3,1,'TileSpacing','tight','Padding','tight')
        Tkt=Tk_diff;
        Tkt = Tkt(Tkt.ebird_n>thr_n&Tkt.name==cities.name(i_c),:);
           
        nexttile([2 1]); hold on; box on; grid on;
        x=(Tkt.takingoff+Tkt.landing);
        y=Tkt.ebird_m;
        f=fit(x,y,"poly1",Weights=Tkt.ebird_n);
        x_ax = floor(min(x)):1:ceil(max(x));
        scatter(x,y,Tkt.ebird_n*2,'.k',"AlphaData",.4)
        p11 = predint(f,x_ax,0.99,'functional','off');
        p = fill([x_ax fliplr(x_ax)]',[p11(:,1);flipud(p11(:,2))],'red', FaceAlpha = .3, EdgeColor = 'none');   
        plot(f)
        xlim([x_ax(1) x_ax(end)])
        xlim([-40 40])
        ylim([-y_lim y_lim])
        xlabel('change on the ground = landing-takeoff'); ylabel('eBird meadian'); title("Corr="+num2str(corr(x,y))); legend off
        
        nexttile; hold on; box on; grid on;
        %x=(Tc.takingoff-Tc.landing)/2;
        x=Tkt.mvt;
        y=abs(Tkt.ebird_m);
        f=fit(x,y,"poly1",Weights=Tkt.ebird_n);
        dx=.01;
        x_ax = (floor(min(x)/dx)*dx):dx:(ceil(max(x)/dx)*dx);
        p11 = predint(f,x_ax,0.99,'functional','off');
        scatter(x,y,Tkt.ebird_n,'.k',"AlphaData",.4)
        p = fill([x_ax fliplr(x_ax)]',[p11(:,1);flipud(p11(:,2))],'red', FaceAlpha = .3, EdgeColor = 'none');   
        plot(f)
        xlim([x_ax(1) x_ax(end)])
        xlabel('max(take-off,landing)'); ylabel('eBird meadian'); title("Corr="+num2str(corr(x,y))); legend off
        ylim([0 y_lim]); yticks(0:2:y_lim)
    end
end
% set(gcf,'renderer','Painters')
% exportgraphics(gcf, "figures_paper/flux_corr_all_spring.eps")

%
figure('position',[0 0 800 350]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
hold on; box on; grid on;
plot([c(:,1) c(:,2)]',[1:height(cities) ; 1:height(cities)],'-k')
scatter(c(:,1),1:height(cities),c(:,3)/300,'s','filled',"MarkerFaceColor",[71 30 15]/255)
scatter(c(:,2),1:height(cities),c(:,3)/300,c(:,2),'filled',"MarkeredgeColor",[71 30 15]/255)
xline(nansum(c(:,1).*c(:,3))/nansum(c(:,3)),"--",num2str(round(1000*nansum(c(:,1).*c(:,3))/nansum(c(:,3)))/10),"color",[71 30 15]/255,"LineWidth",2)
xline(nansum(c(:,2).*c(:,3))/nansum(c(:,3)),"--",num2str(round(1000*nansum(c(:,2).*c(:,3))/nansum(c(:,3)))/10),"color",[210 148 13]/255,"LineWidth",2)
colormap(crameri("batlow"))
set(gca,"ydir","reverse")
xlim([-.2 .2])
ylim([0 height(cities)+1])
yticks(1:height(cities))
yticklabels(cities.name)
xlabel('Correlation')

figure; histogram(Tk_diff.ebird_n)
box on
grid on
xlim([0 100]); xticks([0 50 100])
yticks([0:500:1500])
