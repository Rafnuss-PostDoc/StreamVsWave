

% create a unique id to a pair of chekclit
T2_diff.checklist_id_01 = T2_diff.checklist_id_0 +"_"+ T2_diff.checklist_id_1;

% Load the species level data
load("data/eBird/Tr.mat")
T = vertcat(Tr{:});
Ts = T(:,["species_code","obs_count","checklist_id"]);
Ts = Ts(~isnan(Ts.obs_count),:);


% Filter T to only keep checklist exisit in T2_diff
[id, tmp]= ismember(Ts.checklist_id,T2_diff.checklist_id_0);
Ts0 = Ts(id,:);
Ts0.checklist_id_01 = T2_diff.checklist_id_01(tmp(tmp>0));

[id, tmp]= ismember(Ts.checklist_id,T2_diff.checklist_id_1);
Ts1 = Ts(id,:);
Ts1.checklist_id_01 = T2_diff.checklist_id_01(tmp(tmp>0));

% Join the two daily difference
Ts01 = outerjoin( Ts0, Ts1, "Keys", ["species_code","checklist_id_01"],"MergeKeys",1,"LeftVariables",["obs_count","species_code","checklist_id_01"],"RightVariables",["obs_count"]);

% set nan (no obs) to zero
Ts01.obs_count_Ts0(isnan(Ts01.obs_count_Ts0))=0;
Ts01.obs_count_Ts1(isnan(Ts01.obs_count_Ts1))=0;

% Compute number of individual gained and lost
Ts01.ind_gain = Ts01.obs_count_Ts1-Ts01.obs_count_Ts0;
Ts01.ind_gain(Ts01.ind_gain<0) = 0;
Ts01.ind_lost = Ts01.obs_count_Ts0-Ts01.obs_count_Ts1;
Ts01.ind_lost(Ts01.ind_lost<0) = 0;

% COmpute number of species gained and lost
Ts01.sp_gain = Ts01.obs_count_Ts1>0 & Ts01.obs_count_Ts0==0;
Ts01.sp_lost = Ts01.obs_count_Ts1==0 & Ts01.obs_count_Ts0>0;

% Combine all species
Ts01s = groupsummary(Ts01,"checklist_id_01","sum",["ind_gain", "ind_lost", "sp_gain", "sp_lost"]);

% Join with T2_diff 
T2_diffc = outerjoin( T2_diff, Ts01s, "Keys", "checklist_id_01", "MergeKeys",1);

T2_diffc.sum_sp_gain(isnan(T2_diffc.sum_sp_gain)) = 0;
T2_diffc.sum_sp_lost(isnan(T2_diffc.sum_sp_lost)) = 0;
T2_diffc.sum_ind_gain(isnan(T2_diffc.sum_ind_gain)) = 0;
T2_diffc.sum_ind_lost(isnan(T2_diffc.sum_ind_lost)) = 0;
T2_diffc.GroupCount(isnan(T2_diffc.GroupCount)) = 0;












%% Combine data at the daily scale

% filter 
T2_diffcf = T2_diffc;
T2_diffcf = T2_diffcf(T2_diffcf.cci>-Inf,:);
% T2_diffcf = T2_diffcf(abs(T2_diffcf.ddist)<3,:);
% T2_diffcf = T2_diffcf(abs(T2_diffcf.dhrs)<Inf,:);
% T2_diffcf = T2_diffcf(abs(T2_diffcf.dhours_since_sunrise)<Inf,:);

% group per day (and location)
[G,TID]=findgroups(T2_diffcf(:,["ts","name"]));
% number of species
TID.sp_chg_m=splitapply(@median,T2_diffcf.sum_sp_lost+T2_diffcf.sum_sp_gain,G);
TID.sp_chg_q25=splitapply(@(x) quantile(x,.25),T2_diffcf.sum_sp_lost+T2_diffcf.sum_sp_gain,G);
TID.sp_chg_q75=splitapply(@(x) quantile(x,.75),T2_diffcf.sum_sp_lost+T2_diffcf.sum_sp_gain,G);
%number of individual
TID.ind_chg_m=splitapply(@median,T2_diffcf.sum_ind_lost+T2_diffcf.sum_ind_gain,G);
TID.ind_chg_q25=splitapply(@(x) quantile(x,.25),T2_diffcf.sum_ind_lost+T2_diffcf.sum_ind_gain,G);
TID.ind_chg_q75=splitapply(@(x) quantile(x,.75),T2_diffcf.sum_ind_lost+T2_diffcf.sum_ind_gain,G);

% number of checklist
TID.nb=splitapply(@numel,T2_diffcf.sum_ind_lost+T2_diffcf.sum_ind_gain,G);

% combine with weather radar data
Tk = outerjoin(TID,Trd,MergeKeys=true,RightKeys=["name","obs_dt_day"],LeftKeys=["name","ts"]);
Tk = renamevars(Tk,"ts_obs_dt_day","ts");


% Compute correlation

% filter rain and number per day
thr_nb = 10;
thr_rain = 0; % 5e-6;

c=nan(height(cities),2);
for i_c=1:height(cities)
    Tkt = Tk(Tk.name==cities.name(i_c) & Tk.nb>=thr_nb,:);
    if height(Tkt)>0
        c(i_c,1)=corrW(Tkt.takingoff-Tkt.landing,Tkt.ind_chg_m,Tkt.nb);
        % c(i_c,1)=corrW(Tkt.mvt,Tkt.ind_chg_m,Tkt.nb);
        c(i_c,2) = height(Tkt);
    end
end

figure('position',[0 0 800 350]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
hold on; box on; grid on;
cl = crameri("batlow");
scatter(c(:,1),1:height(cities),c(:,2),'o', "filled", Color=cl(end/2+20,:))
xlim([0 .6])
set(gca,"ydir","reverse")
ylim([0 height(cities)+1])
yticks(1:height(cities))
yticklabels(cities.name)
%legend('Random combinaison','Optimal combinaison','True combinaison',TextColor="w",FontSize=16);
xlabel('Correlation')
%exportgraphics(gcf, "figures/corr_all_spring.png")
%exportgraphics(gcf, "figures_paper/figure_6.png",'BackgroundColor','k')
xline(sum(c(:,1).*c(:,2))/sum(c(:,2)))
nanmean(nansum(c(:,1).*c(:,2))/nansum(c(:,2)))

%%
i_c=3;

figure('position',[0 0 900 750]); tiledlayout(4, 2,'TileSpacing','tight','Padding','tight');

y_wr_max = 200;
for i_y=2021:-1:2017
    Tkt = Tk(year(Tk.ts)==i_y&Tk.name==cities.name(i_c),:);

    if i_y==2021
        nexttile([2,2]);
    else
        nexttile;
    end
    
    hold on; box on; grid on;

    % bar(Tkt.ts,-Tkt.Fd_takingoff(:,i_c)-Fd_landing(:,i_c),1,'FaceColor','w');
    %tmp = cumsum(-Tkt.takingoff-Tkt.landing,'omitnan');
    %tmp = tmp - mean(tmp) + y_wr_max/2;

    bar(Tkt.ts,Tkt.takingoff,1,'FaceColor',clmap(end,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
    bar(Tkt.ts,-Tkt.landing,1,'FaceColor',clmap(1,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
    bar(Tkt.ts,min(Tkt.takingoff,-Tkt.landing),1,FaceAlpha=1,FaceColor=[.85 .85 .85],EdgeColor=[.5 .5 .5])
    % plot(Tkt.ts+.5,tmp,'--k','LineWidth',3)
    % scatter(Tkt.ts(Tkt.rain>5e-6), max(Tkt.takingoff(Tkt.rain>5e-6),-Tkt.landing(Tkt.rain>5e-6)), 50,"k",'filled')
    % scatter(Tkt.ts(Tkt.rain>5e-6), max(Tkt.takingoff(Tkt.rain>5e-6),-Tkt.landing(Tkt.rain>5e-6)), 20,"k",'filled')

    %ylim([0 y_wr_max])
    ylim([0 y_wr_max]);y_lim=ylim;
    yticks(y_lim(1):(y_lim(2)/5):y_lim(2))
    xlim([ datetime(i_y,4,1) datetime(i_y,6,1)])
    text(datetime(i_y,4,1), y_lim(2)," "+num2str(i_y),HorizontalAlignment="left",VerticalAlignment="top",FontSize=24,FontWeight="bold",Color=[.2 .2 .2])
    datetick('x','dd-mmm',"keeplimits")
    if i_y ~= 2021
        xticklabels(''); yticklabels('');
    else
        ylabel("Weather radar nightly change (bird/km^2)")
    end
    set(gca,'XColor','k'); set(gca,'yColor','k');

     yyaxis right

    errorbar(Tkt.ts, Tkt.ind_chg_m, Tkt.ind_chg_m-Tkt.ind_chg_q25, Tkt.ind_chg_q75-Tkt.ind_chg_m,"color",cleBird,"LineStyle","none")
    scatter(Tkt.ts, Tkt.ind_chg_m,Tkt.nb*5,'filled',Color=cleBird)
    ylim([0 40])

end
% exportgraphics(gcf, "figures_paper/figure_5.png",'BackgroundColor','k')









%%

%% Combine with radar

Tq = join(T2_diffc,Trd,RightKeys=["name","obs_dt_day"],LeftKeys=["name","ts"]);


Tqf = Tq;
%Tqf = Tqf(Tqf.cci>-2,:);
% Tqf = Tqf(abs(Tqf.ddist)<1,:);
%Tqf = Tqf(abs(Tqf.dhrs)<Inf,:);
%Tqf = Tqf(abs(Tqf.rain)<=0,:);
%Tqf = Tqf(abs(Tqf.dhours_since_sunrise)<Inf,:);


c=nan(height(cities),2);
for i_c=1:height(cities)
    Tqt = Tqf(Tqf.name==cities.name(i_c),:);
    if height(Tqt)>0
        c(i_c,1)=corr(Tqt.takingoff-Tqt.landing,Tqt.sum_sp_gain+Tqt.sum_sp_lost);
        % corr(Tq.takingoff-Tq.landing,Tq.sum_sp_lost+Tq.sum_sp_gain,"rows","complete");
        %c(i_c,1)=corr(Tkt.takingoff-Tkt.landing,Tkt.ind_chg_m,"rows","complete");
        c(i_c,2) = height(Tqt);
    end
end 

figure('position',[0 0 900 450]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
hold on; box on; grid on;
% scatter(1:height(cities),c(:,1),c(:,2),'o',Color=cl(end/2+20,:), linewidth=2)
plot(1:height(cities),c(:,1),'o')
ylim([-.1 .7])
xlim([0 height(cities)])
xticks(1:height(cities))
xticklabels(cities.name)
%legend('Random combinaison','Optimal combinaison','True combinaison',TextColor="w",FontSize=16);
ylabel('Correlation')
%exportgraphics(gcf, "figures/corr_all_spring.png")
%exportgraphics(gcf, "figures_paper/figure_6.png",'BackgroundColor','k')

yline(nanmean(c(:,1)))
nanmean(c(:,1))

