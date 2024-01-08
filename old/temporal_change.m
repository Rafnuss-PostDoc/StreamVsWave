% Load data generated with download_prepare.m
load("data/eBird/T2r.mat")
load("data/radar_cities")
% load('data/density/inference-trans.mat','g','radar')

% convert to a single table
T2 = vertcat(T2r{:});

% Filter out checklist with high number (active migration count ?)
T2(T2.sum_obs_rcs>3000,:)=[];
% T2(T2.cds_tp>.0005,:)=[];

% Use obs_count of rcs normalized to bird count? 
% Does change much to normalized by rcs but add aditional processing step. 
T2.c = T2.sum_obs_rcs * sum(T2.sum_obs_count) /sum(T2.sum_obs_rcs);
T2.c = T2.sum_obs_count;

% Compute some time variable
T2.weekday = weekday(T2.obs_dt);
T2.obs_dt_day = dateshift(T2.obs_dt,'start','day');
T2.season = (month(T2.obs_dt_day)<7)+1;
T2.doy = day(T2.obs_dt_day,'dayofyear');
T2.year = year(T2.obs_dt_day);

%% Account for effort

effort_var = ["effort_hrs" "effort_distance_km" "cci" "hours_since_sunset" "weekday" "num_observers"];
effort_p_d = [3 4 4 8 0 0];

Mdl = fitglm(T2,"c ~ "+strjoin(effort_var'+"^"+num2str(effort_p_d')," + "));

tmp = T2(1,:);
tmp.effort_hrs=1;
tmp.effort_distance_km=1;
tmp.cci = 1;
tmp.hours_since_sunset=0;
tmp.weekday=3;
tmp.num_observers=1;

% Estimate the number of bird normalized for a standard checklist
T2.cp = Mdl.Residuals.Raw+Mdl.predict(tmp);

% Number can't be below 0
% T2.cp(T2.cp<0)=0;


%% Figure

if false

    figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')
    for i_e=1:numel(effort_var)
        nexttile; hold on; box on; grid on;
        if effort_var{i_e} ~= "num_observers"
            tmp = round(T2.(effort_var{i_e}),1);
        else
            tmp = T2.(effort_var{i_e});
        end
        [G,effort_var_x]=findgroups(tmp);
        tmp = splitapply(@mean,T2.sum_obs_count, G);
        % tmp2 = splitapply(@std,T2.sum_obs_count, G);
        % plot(T2.(effort_var{i_e}),T2.sum_obs_count,'.k');
        plot(effort_var_x,tmp,'.r');
        p=polyfit(T2.(effort_var{i_e}),T2.sum_obs_count,effort_p_d(i_e));
        x=linspace(min(T2.(effort_var{i_e})),max(T2.(effort_var{i_e})),100);
        plot(x,polyval(p,x),'-y','linewidth',2)
        xlabel(effort_var{i_e}); ylabel('sum_obs_rcs')
        axis tight
    end

    %f_poly = @(X,p) [X(:,1).^(1:p(1)-1) X(:,2).^(1:p(2)-1) X(:,3).^(1:p(3)-1) X(:,4).^(1:p(4)-1)];
    %T2fit = [T2.effort_hrs T2.effort_distance_km T2.cci T2.day_of_year];
    %Mdl = fitglm(f_poly(T2fit,effort_p_d),T2.sum_obs_rcs);

    %Mdl = fitrgam(T2,'sum_obs_rcs ~ s(effort_hrs) + effort_distance_km + cci + day_of_year')


    figure; tiledlayout('flow')
    for i_e=1:numel(effort_var)
        nexttile; hold on;
        Mdl.plotPartialDependence(effort_var{i_e})
    end


    figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')
    tmp = sortrows(groupsummary(T2,"name","mean","cci"),"mean_cci","descend");
    for i_name = 1:height(tmp)
        nexttile;
        histogram(T2.cci(T2.name==tmp.name(i_name))); title(tmp.name(i_name)+" M="+tmp.mean_cci(i_name));
        xline(tmp.mean_cci(i_name),'r','LineWidth',2); xlim([-1 2])
    end

    figure; histogram(T2.cp);

    figure;
    for i_name = 1:numel(cityname)
        nexttile; histogram(T2.cp(T2.name==cityname(i_name)))
        xlim([0 200])
    end
    
end


%% Spatial subsample and daily aggregation
if false
    dlatlon=.02; %dlatlon*111
    lon=(q.zone(1)-dlatlon):dlatlon:(q.zone(2)+dlatlon);
    lat=(q.zone(3)-dlatlon):dlatlon:(q.zone(4)+dlatlon);
    [LON,LAT]=meshgrid(lon,lat);


    lon_id = floor((T2.longitude-lon(1)) / dlatlon)+1;
    lat_id = floor((T2.latitude-lat(1)) / dlatlon)+1;
    T2.ll_id = sub2ind([numel(lat), numel(lon)],lat_id,lon_id);

    T2.weekday = weekday(T2.obs_dt);
    % figure; histogram(T2.weekday)

    % Mw = nan(numel(lat), numel(lon),7);
    % nll=numel(lat)*numel(lon);
    % A=groupcounts(T2,{'ll_id','weekday'});
    % Mw(A.ll_id+nll*(A.weekday-1)) = A.GroupCount;

    % figure;tiledlayout('flow','TileSpacing','tight','Padding','tight')
    % for i_w=1:7
    %     nexttile; hold on;
    %     imagesc(lon,lat,Mw(:,:,i_w))
    %     axis tight equal; colorbar;
    % end

    M = nan(numel(lat), numel(lon));
    A=groupcounts(T2,{'ll_id'});
    M(A.ll_id) = A.GroupCount;

    % Select only location with a least thr_M checklists total
    thr_M=1000;
    ll_id_thr = find(M>thr_M);

    figure; hold on;
    % borders('lakes','FaceColor','r')
    plot(maxk(M(:),30),'o-')
    yline(thr_M);  grid on; box on;
    ylabel('Number of checklist'); xlabel('grid cell ordered by number of checklists')

    figure; hold on;
    imagesc(lon,lat,M,'alphadata',M>thr_M)
    plot(T2.longitude,T2.latitude,'.','Color',[.2 .2 .2]);
    plot(LON(ll_id_thr),LAT(ll_id_thr),'.r','MarkerSize',30)
    plot_google_map('MapScale', 1)
    axis tight equal

    T3=T2(ismember(T2.ll_id,ll_id_thr),:);

    % Compute daily summary
    T3.obs_dt_day = dateshift(T3.obs_dt,'start','day');


    % sampling per subgrid
    T4 = groupsummary(T3,{'obs_dt_day','ll_id'},@(x) quantile(x,[.25 .5 .75]),"sum_obs_rcs_res");

    % figure; hold on;
    % plot(T2.obs_dt_day,T2.sum_obs_rcs_res,'.k')
    % plot(T4.obs_dt_day,T4.fun1_sum_obs_rcs_res,'-r')
    % unique_ll=unique(T5.ll_id);
    % for i_ll=1:numel(unique_ll)
    %     T6 = T5(T5.ll_id==unique_ll(i_ll),:);
    %     T6 = sortrows(T6,'obs_dt_day');
    %     plot(T6.obs_dt_day,T6.fun1_sum_obs_rcs_res(:,2))
    % end


    % sampling
    T4 = groupcounts(T3,{'obs_dt_day'});
    T4.fun1_sum_obs_rcs_res=nan(height(T4),3);
    ll_id = unique(T3.ll_id);
    k=5;
    tmp2 = nan(k*numel(ll_id),1);
    for i_du=1:height(T4)
        tmp = T3(T3.obs_dt_day==T4.obs_dt_day(i_du),:);
        for ll_id_i=1:numel(ll_id)
            tmp3=tmp.sum_obs_rcs_res(tmp.ll_id == ll_id(ll_id_i));
            if ~isempty(tmp3)
                tmp2((ll_id_i-1)*k+(1:k)) = datasample(tmp3,k);
            else
                tmp2((ll_id_i-1)*k+(1:k)) = nan;
            end
        end
        T4.fun1_sum_obs_rcs_res(i_du,:) = quantile(tmp2,[.25 .5 .75]);
    end

    % Delete if not more than ... checklists.
    % figure; histogram(T4.GroupCount); box on; grid on; ylabel('Histogram');xlabel('Number of checklists per day')
    T4=T4(T4.GroupCount>3,:);

end

%% Daily aggregation No spatial sampling
    
% compute average number of bird cp, normalize the count with cci. 
% T4 = groupsummary(T2,{'obs_dt_day','name'},@(x) quantile(x,[.25 .5 .75]),"cp");
% T4 = groupsummary(T2,{'obs_dt_day','name'}, {@(x1,x2) sum(x1.*(x2+2))/sum((x2+2)), @(x,y) std(x)/sqrt(length(x))},{"cp", "cci"});
T4 = groupsummary(T2,{'obs_dt_day','name'},{@mean, @var},"cp");

% Sort and compute daily quantity
T4 = sortrows(T4,'obs_dt_day');
T4.season = (month(T4.obs_dt_day)<7)+1;
T4.doy = day(T4.obs_dt_day,'dayofyear');

% Not sure what is that for
% T4doy = groupsummary(T3,{'doy','name'},{@mean, @(x) std(x)/sqrt(length(x))},"cp");


%% Compute timeseries
ty = 1:365;
Fets_mean=nan(numel(ts),numel(cityname));
Fets_var=nan(numel(ts),numel(cityname));
Fets_count=nan(numel(ts),numel(cityname));
Fety = nan(numel(ty),numel(cityname));
for i_name = 1:numel(cityname)
    T4s = T4(T4.name==cityname{i_name},:);
    [Lia,Locb] = ismember(T4s.obs_dt_day,ts);
    %Fets(Locb,:,i_name) = T4s.fun1_sum_obs_count_pred;
    Fets_mean(Locb,i_name) = T4s.fun1_cp;
    Fets_var(Locb,i_name) = T4s.fun2_cp;
    Fets_count(Locb,i_name) = T4s.GroupCount;
    %Fets = Fets-nanmean(Fets);

    T2s = T2(T2.name==cityname{i_name},:);
    [G,ID1]=findgroups(T2s.doy);
    Fety(ID1,i_name) = splitapply(@var,T2s.cp,G);
end


% Smooth estimated variance
TF = isoutlier(Fety,"movmedian",11,"ThresholdFactor",5);
Fety2 = Fety; Fety2(TF)=nan;
Fety3 = reshape(smooth(Fety2,11),size(Fety2));
Fety3(all(isnan(Fety),2),:)=nan;

figure; tiledlayout(3,1,'TileSpacing','none','Padding','none')
nexttile; plot(Fety); axis([50 350 0 1400])
nexttile; plot(Fety2); axis([50 350 0 1400])
nexttile; plot(Fety3); axis([50 350 0 1400])

% Variance of estimation with law of large number
tmp = day(ts,'dayofyear'); tmp(tmp==366)=365;
Fets_mean_var = Fety3(tmp,:)./Fets_count;

%% Filter to keep only if we have enough checklist
n_thr=10;

figure;
imagesc(datenum(ts),1:numel(cityname),(Fets_count<n_thr)')
yticks(1:numel(cityname)); yticklabels(cityname); datetick('x'); % caxis([0 40])

Fets_mean(Fets_count<n_thr)=nan;
Fets_var(Fets_count<n_thr)=nan;

if false
    % T-test
    tmp = (Fets_count-1).*Fets_var;
    dof = Fets_count(1:end-1,:) + Fets_count(2:end,:) - 2 ;
    tmp2 = (tmp(1:end-1,:)+tmp(2:end,:) ) ./ dof;
    tmp3 = sqrt(tmp2) .* sqrt(1./Fets_count(1:end-1,:) + 1./Fets_count(2:end,:));
    t = ( Fets_mean(1:end-1,:)-Fets_mean(2:end,:) ) ./ tmp3;
    
    p = tcdf(t,dof);
    
    figure;
    imagesc(datenum(ts(1:end-1)),1:numel(cityname),p')
    yticks(1:numel(cityname)); yticklabels(cityname); datetick('x'); % caxis([0 40])
end


%% Average seasonal pattern
figure; tiledlayout('flow','TileSpacing','tight','Padding','tight'); set(gcf, 'color', 'k');
for i_name=1:numel(cityname)
    nexttile;hold on; box on; grid on; hold on;
    title(cityname{i_name},'Color','w')
    set(gca,'color', 'k');
    set(gca,'XColor','w');
    set(gca,'YColor','w');
    plot(Fety(:,i_name))
    % plot(Fety(:,i_name)+Fety3(:,i_name),'--')
    % plot(Fety(:,i_name)-Fety3(:,i_name),'--')
    ylabel('eBird')
    yyaxis right; box on
    tmp = -Fd_takingoff(:,i_name)-Fd_landing(:,i_name);
    plot(median(reshape(tmp(1:(365*11)),365,[]),2))
    ylabel('Weather Radar')
end




%% Radar 
col2 = [0.4660 0.6740 0.188];
col3 = [0.9290    0.6940    0.1250];
clmap = crameri('berlin');
ts_season = (month(ts)<7)+1;

for i_y=2010:2020
    figure('position',[0 0 1650 1200]); tiledlayout('flow','TileSpacing','tight','Padding','tight'); set(gcf, 'color', 'k');
    for i_name=1:numel(cityname)
        nexttile;hold on; box on; grid on; title(cityname{i_name},'color','w')
        set(gca, 'color', 'k');
        set(gca,'XColor','w');
        set(gca,'yColor','w');
        bar(ts,-Fd_takingoff(:,i_name),1,'FaceColor',clmap(end,:))
        bar(ts,-Fd_landing(:,i_name),1,'FaceColor',clmap(1,:))
        t=bar(ts,-Fd_takingoff(:,i_name)-Fd_landing(:,i_name),1,'FaceColor','w');
        plot(ts,Fd_diff(:,i_name),'linewidth',2,'color','w')
        %ylabel('Bird density (bird/km^2)')
        xlim([datetime(i_y,4,1) datetime(i_y,11,15)])
        ylim([-1 1]*800)
        xticklabels(''); yticklabels('');

    end
    % exportgraphics(gcf, "figures/presentation/comulative_all_"+num2str(i_y)+".png",'BackgroundColor','k')
end


%% Radar vs eBird

% i_name=2;

figure('position',[0 0 1600 900]);
    % tiledlayout(6,10,'TileSpacing','tight','Padding','tight')
    % tiledlayout(3,2,'TileSpacing','tight','Padding','tight')
    tiledlayout(4,3,'TileSpacing','tight','Padding','tight')
for i_y=2020% 2015:2020
 

for i_name=1:numel(cityname)
    T3s = T3(T3.name==cityname{i_name},:);
    
    for i_s=2%1:2

        y_lim_r= [-220 220];% [min(Fd_diff(ts_season==i_s,i_name)) max(Fd_diff(ts_season==i_s,i_name))];
        y_lim_e = [-20 20];%[100 200]; % [min(Fets(ts_season==i_s,2,i_name)) max(Fets(ts_season==i_s,2,i_name))];

       
            nexttile;hold on;

            id1=year(T3s.obs_dt_day)==i_y & T3s.season==i_s;
            plot(T3s.obs_dt_day(id1),T3s.cp(id1),'.','markersize',10,'color',[.3 .3 .3])

            id2 = find(year(ts)==i_y & ts_season==i_s);
            plot(ts(id2),Fets_mean(id2,i_name),'-','linewidth',2,'color',col2)
            plot(ts(id2),Fets_mean(id2,i_name)+Fets_mean_var(id2,i_name),'--','color',col2)
            plot(ts(id2),Fets_mean(id2,i_name)-Fets_mean_var(id2,i_name),'--','color',col2)
            scatter(ts(id2),Fets_mean(id2,i_name),Fets_count(id2,i_name)/3,col2,'filled')

            ylim(nanmean(Fets(id2,1,i_name))+y_lim_e)


            yticklabels('')
            yyaxis right
            plot(ts(id2),Fd_diff(id2,i_name)-Fd_diff(id2(1),i_name),'-','linewidth',2,'color','w')
            tmp = mean(ylim);
            ylim(tmp+y_lim_r)
            if i_s==2
                xlim([datetime("1-April-"+num2str(i_y)) datetime("1-June-"+num2str(i_y))])
                title(cityname{i_name},"color","w");%+" | Spring")
                %title(cityname{i_name}+" | "+i_y,"color","w")
            else
                xlim([datetime("1-Aug-"+num2str(i_y)) datetime("15-Nov-"+num2str(i_y))])
                % title(cityname{i_name}+" | Autumn")
            end


            yticklabels('')
            xticklabels('')
            box on; grid off;
            set(gca,'Color','k')
            ax = gca;
            ax.YAxis(1).Color = 'w';
            ax.YAxis(2).Color = 'w';
            ax.XAxis.Color = 'w';
        end
    end
    % exportgraphics(gcf, "figures/presentation/ebird_"+cityname(i_name)+"_spring.png",'BackgroundColor','k')

end


%%

figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')

r22=nan(numel(cityname),6);

for i_name=1:numel(cityname)
    for i_y=2015:2020
        id2 = find(year(ts)'==i_y & ts_season'==i_s & ~isnan(Fets(:,1,i_name)) & ~isnan(Fd_diff(:,i_name)));
        r2=corrcoef(Fets(id2,1,i_name),Fd_diff(id2,i_name));
        r22(i_name,i_y-2014)=r2(1,2);
    end
end


for a=1:2

    nexttile; hold on;
    x=Fsall(:,i_y,1);
    y=Fsall(:,i_y,3);

    plot(x,y,'.k','markersize',10);
    for i_y2=1:11
        y2=Fsall(:,i_y2,3);
        plot(x,y2,'.','Color',[.7 .7 .7]);
        r2=corrcoef(y2,x,'Row','complete');
        r22(i_y2)=r2(1,2);
    end
    l = lsline;
    plot(x,y,'.k','markersize',10);
    r2=corrcoef(y(:),x(:),'Row','complete');
    legend(l(end:1-end), "R^2="+num2str(r2(1,2)), "R^2="+num2str(nanmean(r22))+" ("+num2str(nanstd(r22))+")")
    axis equal tight square; box on; title(num2str(i_y))
end


%%
lm1 = cell(11,11);
r2 = nan(11,11);
rmse = nan(11,11);

% lm2 = fitlm([reshape(Fsall(1:end-5,:,3),[],1) repmat(t(1:end-5),11,1)],reshape(Fsall(1:end-5,:,1),[],1),'Weights',reshape(Fsall(1:end-5,:,5),[],1));

figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')

t=(1:size(Fsall,1))';
for i_y = 1:11
    nexttile; hold on
    x=Fsall(:,i_y,1);
    % plot(x)

    e = Fsall(:,:,3)-x;

    T =repmat(t,11,1);
    lm3 = fitlm(T,e(:));

    e = e-reshape(lm3.predict(T),size(e));

    rmse(i_y,:) = nanmean(e.^2);

    plot(e)
    plot(e(:,i_y),'r','linewidth',3)
    yline(0,'--k')
    %for i_y2 = 1:11
    %y=Fsall(:,i_y2,3);
    %w=Fsall(:,i_y2,5);%1./abs(Fsall(:,i_y2,4)-Fsall(:,i_y2,1));

    %         id=~isnan(x)&~isnan(y);
    %         x=x(id);
    %         y=y(id);

    % lm1{i_y,i_y2} = fitlm([y t],x,'Weights',w);
    % plot(lm1{i_y,i_y2}.predict(y),'.')
    %plot(lm1{i_y,i_y2}.predict([y t]),'.')
    %r2(i_y,i_y2) = lm1{i_y,i_y2}.Rsquared.Adjusted;
    %rmse(i_y,i_y2) = lm1{i_y,i_y2}.RMSE;

    % xpred = lm2.predict([y t]);

    % e=x-lm2.predict([y t]);

    %         e=x-y; e-nan
    %
    %         plot(e,'.')
    %         rmse(i_y,i_y2) = nanmean(e(60:150).^2);
    %         axis tight;
    xlim([60 150]); ylim([-200 200])
    %end
end
% figure; imagesc(rmse./mean(rmse,2))




%%
figure; tiledlayout('flow','Padding','none','TileSpacing','compact')
nexttile; x=Fets(1:end,2); y=mean(abs(Fts),2);
plot(x, y,'.k'); l =lsline; l.Color='r'; r2=corrcoef(y(:),x(:),'Row','complete'); legend("R^2="+num2str(r2(1,2)))
xlabel('bird on the ground'); ylabel('mean abs(land. dep.) during previous night')

nexttile; x=diff(Fets(:,2)); y=mean(abs(Fts(2:end,:)),2);
plot(x, y,'.k'); l =lsline; l.Color='r'; r2=corrcoef(y(:),x(:),'Row','complete'); legend("R^2="+num2str(r2(1,2)))
xlabel('Diff bird on the ground from previous day'); ylabel('mean land. dep. during previous night')

nexttile; x=Fets(:,2); y=sum(Fts,2);
plot(x, y,'.k'); l =lsline; l.Color='r'; r2=corrcoef(y(:),x(:),'Row','complete'); legend("R^2="+num2str(r2(1,2)))
xlabel(' bird on the ground '); ylabel('diff land. dep. during previous night')

nexttile; x=diff(Fets(:,2)); y=sum(Fts(2:end,:),2);
plot(x, y,'.k'); l =lsline; l.Color='r'; r2=corrcoef(y(:),x(:),'Row','complete'); legend("R^2="+num2str(r2(1,2)))
xlabel(' Diff bird on the ground '); ylabel('diff land. dep. during previous night')

nexttile; x=abs(Fets); y=(sum(Fts,2));
plot(x, y,'.k'); l =lsline; l.Color='r'; r2=corrcoef(y(:),x(:),'Row','complete'); legend("R^2="+num2str(r2(1,2)))
xlabel(' abs( bird on the ground) '); ylabel('abs(diff land. dep.) during previous night')

