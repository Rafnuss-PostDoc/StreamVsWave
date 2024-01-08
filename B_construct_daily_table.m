% Load data generated with download_prepare.m
load("data/eBird/T2r.mat")
load("data/radar_cities")
load("data/rain/cities_rain","rain")
clmap = crameri('berlin');
cleBird = [36 128 52]/255;
cleBird2 = [240 246 250]/255;
% load('data/density/inference-trans.mat','g','radar')

% convert to a single table
T2 = vertcat(T2r{:});

% Compute some time variable
T2.weekend(:) = "weekday";
T2.weekend(weekday(T2.obs_dt)==1 | weekday(T2.obs_dt)==7)="weekend";
T2.obs_dt_day = dateshift(T2.obs_dt,'start','day');
% figure; histogram(month(T2.obs_dt))
T2.doy = day(T2.obs_dt_day,'dayofyear');

T2.year = year(T2.obs_dt_day);

% Use obs_count of rcs normalized to bird count?
% Does change much to normalized by rcs but add aditional processing step.
T2.c = T2.sum_obs_rcs * sum(T2.sum_obs_count) /sum(T2.sum_obs_rcs);
T2.c = T2.sum_obs_count;


%% Filter

% filtering time of day
id = T2.doy > 150;
T2(id,:)=[];

% Filter out checklist with high number (active migration count ?)
% T2(T2.sum_obs_rcs>3000,:)=[];
%T2(T2.sum_obs_count>500,:)=[]; % slightly improve but not much

% filter rain
% mean(T2.cds_tp>0) 1.6%
id = T2.cds_tp>0;
disp(num2str(round(mean(id)*100,2)) + "% of dataset flagged as raining")
% T2(id,:)=[]; remove
T2.precipitation = id;

% filtering time of day
id = T2.ss~=0;
disp(num2str(round(mean(id)*100,2)) + "% of dataset remove because during night")
T2(id,:)=[];

% Protocol
id = T2.protocol_id=="P22" | T2.protocol_id=="P21";
disp(num2str(round(mean(id)*100,2)) + "% of dataset is stationary or traveling")
T2=T2(id,:);

% Filtering CCI
id = T2.cci<-2.5;
disp(num2str(round(mean(id)*100,2)) + "% of dataset has cci<-2.5. remove")
T2(id,:)=[];
id = T2.cci>2.5;
disp(num2str(round(mean(id)*100,2)) + "% of dataset has cci>2.5. Threashold to 2.5")
T2.cci(id)=2.5;

% num of observator
id = T2.num_observers>20;
disp(num2str(round(mean(id)*100,2)) + "% of dataset has more than 20 observers. Threashold to 20")
T2.num_observers(id)=20;

% On the exact hours
T2.effort_hrs_on_hour = round(T2.effort_hrs)==T2.effort_hrs;
disp(num2str(round(mean(T2.effort_hrs_on_hour)*100,2)) + "% of dataset is on the hour (duration is 1,2,3... hours)")

% figure; nexttile; histogram(T2.c); nexttile; histogram(log(T2.c))

%% Account for effort

effort_var = ["effort_hrs" "effort_distance_km" "cci" "hours_since_sunrise" "num_observers" "cds_i10fg" "cds_msl" "cds_lcc" "precipitation" "effort_hrs_on_hour" "protocol_id" "weekend"];
effort_p_d = [3              3                    3                4          3                 3        2            1                 1               1            1              1];


if false
    figure('position',[0 0 800 700]); tiledlayout(4,4,'TileSpacing','tight','Padding','tight')
    cm=colormap;
    cm_f = @(a,b,v) cm(min(max(1,round((v-a)./(b-a)*height(cm))),height(cm)),:);
    for i_e=1:numel(effort_var)
        nexttile; hold on; box on; grid on;
        if isstring(T2.(effort_var{i_e}))
            tmp = categorical(T2.(effort_var{i_e}));
        else
            rd = -round(log10(range(T2.(effort_var{i_e}))/100));
            if rd==Inf, rd=1;end;
            tmp = round(double(T2.(effort_var{i_e})),rd);
        end
        [G,effort_var_x] = findgroups(tmp);
        tmp2 = splitapply(@mean,T2.c, G);
        n = splitapply(@numel,T2.c, G);
        
        if i_e<9
            scatter(effort_var_x,tmp2,100,n,'filled');%colorbar;
            p=fit(tmp,T2.c,"poly"+effort_p_d(i_e));
            x=linspace(min(T2.(effort_var{i_e})),max(T2.(effort_var{i_e})),100);
            plot(x,p(x),'-r','linewidth',2)
        else
            b=bar(categorical(effort_var_x),tmp2,'FaceColor','flat');
            b.CData = cm_f(0,sum(n),n);
        end
        xlabel(effort_var{i_e},'Interpreter','none'); %ylabel('sum_obs_rcs')
        axis tight
        ylim([0 80]);
    end

    nexttile([1 4]); hold on; box on; grid on;
    [G,effort_var_x] = findgroups(T2.name);
    tmp = splitapply(@mean,T2.c, G);
    n = splitapply(@numel,T2.c, G);
    scatter(categorical(effort_var_x, cities.name),tmp, n/100, n,'filled')
    % b=bar(categorical(effort_var_x, cities.name),tmp,'FaceColor','flat');%colorbar;
    % b.CData = cm_f(min(n),max(n),n);
    xlabel("City"); ylabel('sum_obs_rcs')
    h=text(b.XEndPoints,b.YEndPoints+3,effort_var_x,'HorizontalAlignment','left','VerticalAlignment','middle')
    set(h,'Rotation',90); xticklabels([])
    axis tight
    ylim([0 80]);
end


% Mdl =  fitlme(T2,"c ~ "+strjoin(effort_var'+"^"+num2str(effort_p_d')," + ")+" + effort_hrs_on_hour + (1|protocol_id) + (1|name) + (1|weekend)");

T2t = T2(:,["c" effort_var "name"]);
T2t.cds_msl = T2t.cds_msl/1e5;

Mdl = fitglme(T2t,"c ~ "+strjoin(effort_var'+"^"+num2str(effort_p_d')," + ")+" + effort_hrs_on_hour^1 + protocol_id^1 + weekend^1 + (1|name)", Distribution="poisson");

% Mdl2 = fitrensemble(T2t,"c"); % NPrint=1 slow

tmp = T2t;
tmp.effort_hrs(:)=1; % mean 1.2
tmp.effort_hrs_on_hour(:)=0; 
tmp.effort_distance_km(:)=1; % mean=0
tmp.cci(:) = 1;
tmp.hours_since_sunrise(:)=2; % mean=6
tmp.weekend(:)="weekday";
tmp.num_observers(:)=1; % mean=1.6
tmp.protocol_id(:)="P22"; % traveling
tmp.cds_i10fg(:) = 5; % mean(T2t.cds_i10fg)=8
tmp.cds_msl(:) = mean(T2t.cds_msl); % 1.0158
tmp.cds_lcc(:) = 0; % mean 0.25
tmp.precipication(:) = 0;

% Mdl.predict(tmp)

% Estimate the number of bird normalized for a standard checklist
T2.cp = Mdl.Residuals.Raw+Mdl.predict(tmp);

% T2.cp = T2t.c-Mdl2.predict(T2t)+Mdl2.predict(tmp);

%% Figure Effort

if false 
    
    figure('position',[0 0 1200 350]); tiledlayout(1,numel(effort_var)+2,'TileSpacing','tight','Padding','tight')
    for i_e=1:numel(effort_var)
        nexttile; hold on;
        Mdl.plotPartialDependence(effort_var{i_e}, NumObservationsToSample=0)
        ylim([0 100]);
        if i_e==1, ylabel('Total number of bird seen per checklists'); else, ylabel(""); yticklabels([]), end
        title("")
        box on; grid on
    end
    nexttile; hold on;
    for i_c=1:height(cities)
        Mdl{i_c}.plotPartialDependence("weekend", NumObservationsToSample=0)
    end
    title(""); box on; grid on; ylim([0 100]);
    nexttile; hold on;
    for i_c=1:height(cities)
        Mdl{i_c}.plotPartialDependence("city", NumObservationsToSample=0)
    end
    title(""); box on; grid on; ylim([0 100]);
    
    
    figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')
    tmp = sortrows(groupsummary(T2,"name","mean","cci"),"mean_cci","descend");
    for i_c = 1:height(tmp)
        nexttile;
        histogram(T2.cci(T2.name==tmp.name(i_c))); title(tmp.name(i_c)+" M="+tmp.mean_cci(i_c));
        xline(tmp.mean_cci(i_c),'r','LineWidth',2); xlim([-1 2])
    end
    
    figure; histogram(T2.cp);
    
    figure;
    for i_c = 1:height(cities)
        nexttile; histogram(log(T2.c(T2.name==cities.name(i_c))))
        xlim([0 7])
    end
end

%% Combine all data in a table

Trd = table(...
    repelem(cities.name,numel(ts),1),...
    repmat(ts',height(cities),1),...
    Fd_takingoff(:),...
    Fd_landing(:),...
    Fd_mvt(:),...
    rain(:),...
    VariableNames=["name","obs_dt_day","takingoff","landing","mvt","rain"]);

[G,TID]=findgroups(T2(:,["obs_dt_day","name"]));
TID.ebird_n=splitapply(@numel,T2.cp,G);
TID.ebird_q10=splitapply(@(x) quantile(x,.10),T2.cp,G);
TID.ebird_q25=splitapply(@(x) quantile(x,.25),T2.cp,G);
TID.ebird_q50=splitapply(@(x) quantile(x,.5),T2.cp,G);
TID.ebird_q75=splitapply(@(x) quantile(x,.75),T2.cp,G);
TID.ebird_q90=splitapply(@(x) quantile(x,.9),T2.cp,G);
TID.ebird_d=splitapply(@(x) {x} ,T2.cp,G);

Tk = outerjoin(TID,Trd,MergeKeys=true,RightKeys=["name","obs_dt_day"],LeftKeys=["name","obs_dt_day"]);
Tk = renamevars(Tk,"obs_dt_day","ts");





