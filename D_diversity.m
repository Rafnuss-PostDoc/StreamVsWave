

%
% run("B_construct_daily_table.m")

%% Load the species level data
load("data/eBird/Tr.mat")
T = vertcat(Tr{:});
Ts = T(:,["species_code","obs_count","checklist_id", "obs_dt_day", "name"]);

% Filter by checklist belonging to T2
[id, tmp]= ismember(Ts.checklist_id,T2.checklist_id);
Ts = Ts(id,:);
Ts = renamevars(Ts,"obs_dt_day","ts");

% Find counts of all species per day and location
Tsdns = groupsummary(Ts,["ts", "name","species_code"],'sum','obs_count');
Tsdns.GroupCount(Tsdns.sum_obs_count==0)=0;

Tsdns_sp = unstack(Tsdns(:,[1:3 5]),"sum_obs_count","species_code");
Tsp = Tsdns_sp(:,1:2);
Tsp.nb_sp = sum(Tsdns_sp{:,3:end}>0,2,'omitnan');
Tsp.nb_ind = sum(Tsdns_sp{:,3:end},2,'omitnan');

% Join with the other variable
Tc = innerjoin(Tk, Tsp, keys=["ts","name"],LeftVariables=["ts","rain","name","ebird_n","landing","takingoff","mvt"]);






%% Change in number of bird species
% not very intresting... 

i_c=3;
figure('position',[0 0 900 750]); tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
y_eb_max = 40+nanmean(Tc.nb_sp(Tc.name==cities.name(i_c)));
y_eb_min = -40+nanmean(Tc.nb_sp(Tc.name==cities.name(i_c)));
y_wr_max = (y_eb_max-y_eb_min)*5;
for i_y=2021:-1:2017
    Tct = Tc(year(Tc.ts)==i_y&Tc.name==cities.name(i_c)&Tc.ebird_n>0,:);
    if i_y==2021
        nexttile([2,2]);
    else
        nexttile;
    end
    hold on; box on; grid on;

    plot(Tct.ts+.5, Tct.nb_sp, Color=cleBird, linewidth=2)
   
    ylim([y_eb_min y_eb_max])

    if i_y ~= 2021
        xticklabels(''); yticklabels('');
    else
        ylabel("NUmber of species")
    end
    yyaxis right

    % bar(Tct.ts,-Tct.Fd_takingoff(:,i_c)-Fd_landing(:,i_c),1,'FaceColor','w');
    tmp = cumsum(-Tct.takingoff-Tct.landing,'omitnan');
    tmp = tmp - mean(tmp) + y_wr_max/2;

    bar(Tct.ts,Tct.takingoff,1,'FaceColor',clmap(end,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
    bar(Tct.ts,-Tct.landing,1,'FaceColor',clmap(1,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
    bar(Tct.ts,min(Tct.takingoff,-Tct.landing),1,FaceAlpha=1,FaceColor=[.85 .85 .85],EdgeColor=[.5 .5 .5])
    plot(Tct.ts+.5,tmp,'--k','LineWidth',3)
    % scatter(Tct.ts(Tct.rain>5e-6), max(Tct.takingoff(Tct.rain>5e-6),-Tct.landing(Tct.rain>5e-6)), 50,"k",'filled')
    % scatter(Tct.ts(Tct.rain>5e-6), max(Tct.takingoff(Tct.rain>5e-6),-Tct.landing(Tct.rain>5e-6)), 20,"k",'filled')

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
end


%% 
thr_n=50;
thr_rain = 5e-6;
c = nan(height(cities),3);

for i_c=1:height(cities)
    Tct=Tc;
    Tct = Tct(Tct.ebird_n>thr_n&Tct.name==cities.name(i_c),:);
    Tct(Tct.rain>=thr_rain,:)=[];
    Tctoff=Tct;
    Tctoff.ts = Tctoff.ts+1;
    Tctd = innerjoin(Tct,Tctoff,keys=["ts","name"],RightVariables=["nb_sp","ebird_n"]);
    Tctd = renamevars(Tctd,["nb_sp_Tct","nb_sp_Tctoff"],["nb_sp_0","nb_sp_1"]);
    Tctd.nb_sp_01 = Tctd.nb_sp_1-Tctd.nb_sp_0;
    Tctd.ebird_n = min([Tctd.ebird_n_Tct Tctd.ebird_n_Tctoff],[],2);

    if sum(Tctd.ebird_n)>0
        c(i_c,3) = sum(Tctd.ebird_n);
        c(i_c,1) = corrW(Tctd.mvt,abs(Tctd.nb_sp_01),Tctd.ebird_n);
        c(i_c,2) = corrW(Tctd.takingoff+Tctd.landing,Tctd.nb_sp_01,Tctd.ebird_n);
    end
end

figure('position',[0 0 800 350]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
hold on; box on; grid on;
xline(nansum(c(:,1).*c(:,3))/nansum(c(:,3)),"--r",num2str(nansum(c(:,1).*c(:,3))/nansum(c(:,3))))
xline(nansum(c(:,2).*c(:,3))/nansum(c(:,3)),"--k",num2str(nansum(c(:,2).*c(:,3))/nansum(c(:,3))))
scatter(c(:,1),1:height(cities),c(:,3)/300,"r",'filled')
scatter(c(:,2),1:height(cities),c(:,3)/300,'k','filled')

colormap(crameri("batlow"))
set(gca,"ydir","reverse")
% xlim([0 .5])
axis tight
ylim([0 height(cities)+1])
yticks(1:height(cities))
yticklabels(cities.name)
legend('Mvt vs abs(ebird change)','WR change on the ground vs ebird change',"n_thr="+num2str(thr_n),location="northoutside",Orientation="horizontal");
xlabel('Correlation')





















%% Frequency of observation

Tsdns_obs = unstack(Tsdns(:,1:4),"GroupCount","species_code");

% Join with the other variable
Tc = innerjoin(Tk, Tsdns_obs, keys=["ts","name"],LeftVariables=["ts","rain","name","ebird_n","landing","takingoff","mvt"]);

% frequency of observation (number of observation / number of checklist)
M=Tc{:,8:end}./Tc.ebird_n;
M(isnan(M))=0;

Tc = Tc(:,1:7);
Tc.beta(:) = nan;
Tc.beta_a(:) = nan;
Tc.beta_b(:) = nan;

for i_c=1:height(cities)
    for i_y=2021:-1:2017
        id = find((year(Tc.ts)==i_y) & cities.name(i_c)==Tc.name);
        assert(all(hours(diff(Tc.ts(id)))==24))
        
% Keep this because it is wrong, but a nice try! 
% you can use the ratio of report as a measure of species presence when computing the difference
% The expected change in species is either 0-1 or 0-1 (both 1-1 and 0-0
% would not count for a change in species). 

%         % compute the change in ratio of report of each species and sum the
%         % abs value of all species for a given days
%         U = sum(abs(diff(M(id,:))),2);
% 
%         % sort of richness between consecutive, normalize the change in
%         % diversity between the two days. 
%         % not sure we should be using this. I think it's a matter of if
%         % later in the year, more species in fewer number occurs or not
%         % really. 
%         D = sum(max(cat(3,M(id(1:end-1),:),M(id(2:end),:)),[],3),2);
% 
%         % assigne between consecutive day at the mid-nighte between the two
%         % days. (so difference day 1 and day 2 is assign at midnight day
%         % 2). 
%         Tc.beta(id(2:end))=U./D;


        a=M(id(1:end-1),:);
        b=M(id(2:end),:);

        % we compute the exptectation that we have a species in a and not
        % in b + the species in b but not in a (ie., a change in the number
        % of species in a checklists). 
        % This is equivalent to drawing randomly two checklists and
        % computing the sum of species change. 
        U = sum(a.*(1-b)+(1-a).*b,2);
        % for D we compute the probability of in a or b or a and b (so 1-
        % neither a nor b). 
        D = sum(1-((1-a).*(1-b)),2);

        Tc.beta(id(2:end)) = U; % ./D;

        Tc.beta_a(id(2:end)) = sum(a.*(1-b),2);
        Tc.beta_b(id(2:end)) = sum((1-a).*b,2);
    end
end



%% 

i_c=3;
figure('position',[0 0 900 750]); tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
y_eb_max = 10+nanmean(Tc.beta(Tc.name==cities.name(i_c)));
y_eb_min = -10+nanmean(Tc.beta(Tc.name==cities.name(i_c)));
y_wr_max = (y_eb_max-y_eb_min)*10;
for i_y=2021:-1:2017
    Tct = Tc(year(Tc.ts)==i_y&Tc.name==cities.name(i_c)&Tc.ebird_n>0,:);
    if i_y==2021
        nexttile([2,2]);
    else
        nexttile;
    end
    hold on; box on; grid on;

    plot(Tct.ts, Tct.beta, Color=cleBird, linewidth=2)
   
    ylim([y_eb_min y_eb_max])

    if i_y ~= 2021
        xticklabels(''); yticklabels('');
    else
        ylabel("NUmber of species")
    end
    yyaxis right

    % bar(Tct.ts,-Tct.Fd_takingoff(:,i_c)-Fd_landing(:,i_c),1,'FaceColor','w');
    tmp = cumsum(-Tct.takingoff-Tct.landing,'omitnan');
    tmp = tmp - mean(tmp) + y_wr_max/2;

    bar(Tct.ts,Tct.takingoff,1,'FaceColor',clmap(end,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
    bar(Tct.ts,-Tct.landing,1,'FaceColor',clmap(1,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
    bar(Tct.ts,min(Tct.takingoff,-Tct.landing),1,FaceAlpha=1,FaceColor=[.85 .85 .85],EdgeColor=[.5 .5 .5])
    plot(Tct.ts+.5,tmp,'--k','LineWidth',3)
    % scatter(Tct.ts(Tct.rain>5e-6), max(Tct.takingoff(Tct.rain>5e-6),-Tct.landing(Tct.rain>5e-6)), 50,"k",'filled')
    % scatter(Tct.ts(Tct.rain>5e-6), max(Tct.takingoff(Tct.rain>5e-6),-Tct.landing(Tct.rain>5e-6)), 20,"k",'filled')

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
end


%% 
thr_n=0;
thr_rain = 5e-6;
c = nan(height(cities),3);

for i_c=1:height(cities)
    Tct=Tc(Tc.ebird_n>thr_n&Tc.name==cities.name(i_c),:);
    Tct(Tct.rain>=thr_rain,:)=[];
    
    if sum(Tct.ebird_n)>0
        c(i_c,3) = sum(Tct.ebird_n);
        c(i_c,1) = corrW(Tct.mvt,Tct.beta,Tct.ebird_n);
        c(i_c,2) = corrW(Tct.takingoff+abs(Tct.landing),Tct.beta,Tct.ebird_n);
    end
end

figure('position',[0 0 800 350]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
hold on; box on; grid on;
plot([c(:,1) c(:,2)]',[1:height(cities) ; 1:height(cities)],'-k')
scatter(c(:,1),1:height(cities),c(:,3)/300,'s','filled',"MarkerFaceColor",[71 30 15]/255)
scatter(c(:,2),1:height(cities),c(:,3)/300,c(:,2),'filled',"MarkeredgeColor",[71 30 15]/255)
xline(nansum(c(:,1).*c(:,3))/nansum(c(:,3)),"--",num2str(round(1000*nansum(c(:,1).*c(:,3))/nansum(c(:,3)))/10),"color",[71 30 15]/255,"LineWidth",2)
xline(nansum(c(:,2).*c(:,3))/nansum(c(:,3)),"--",num2str(round(1000*nansum(c(:,2).*c(:,3))/nansum(c(:,3)))/10),"color",[210 148 13]/255,"LineWidth",2)
colormap(crameri("batlow"))
set(gca,"ydir","reverse")
xlim([0 1])
ylim([0 height(cities)+1])
yticks(1:height(cities))
yticklabels(cities.name)
xlabel('Correlation')
