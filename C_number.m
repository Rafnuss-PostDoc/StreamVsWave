
% run B_construct_daily_table
path_bmmus="/Users/rafnuss/Library/CloudStorage/Box-Box/BMM-US";
load(path_bmmus+'/data/density/inference-trans.mat','g','radar')

y=2010:2021;
Tk=Tk(Tk.ebird_n>0,:);
%% Daily correlation

% remove what appear as outliar
thr_rain = 5e-6; % rain
thr_n = 50; % no enough checklists
Tkfit = Tk(Tk.ebird_q75<120 & Tk.rain<=thr_rain & Tk.ebird_n>=thr_n,:);

c = nan(height(cities),3);
fit_chg_ebird = nan(height(cities),3);
fit_chg_ebird_y = nan(height(cities),numel(y),2);

for i_c=1:height(cities)
    Tkt = Tkfit(Tkfit.name==cities.name(i_c),:);
    Tktoff = Tkt;
    Tktoff.ts = Tktoff.ts+1;
    Tc = innerjoin(Tkt,Tktoff,keys=["ts","name"],RightVariables=["ebird_q75","ebird_n"]);
    Tc = renamevars(Tc,["ebird_q75_Tkt","ebird_q75_Tktoff"],["ebird_q75_0","ebird_q75_1"]);
    Tc.ebird_q75_01 = Tc.ebird_q75_1-Tc.ebird_q75_0;
    Tc.ebird_n = sum([Tc.ebird_n_Tkt Tc.ebird_n_Tktoff],2);

    if sum(Tc.ebird_n)>0
        c(i_c,3) = sum(Tc.ebird_n);
        c(i_c,1) = corrW(Tc.mvt,abs(Tc.ebird_q75_01),(Tc.ebird_n));
        c(i_c,2) = corrW(Tc.takingoff+Tc.landing,Tc.ebird_q75_01,(Tc.ebird_n));
        % c(i_c,1) = corr(Tc.mvt,abs(Tc.ebird_q75_01));
        % c(i_c,2) = corr(Tc.takingoff+Tc.landing,Tc.ebird_q75_01);
        c(i_c,3) = sum(Tc.ebird_n);
        tmp = fit(Tc.takingoff+Tc.landing, Tc.ebird_q75_01, fittype('a*x'));
        fit_chg_ebird(i_c,:) = 1./[tmp.a; confint(tmp,.90)]';
        for i_y=1:numel(y)
            Tc2 = Tc(year(Tc.ts)==y(i_y),:);
            if height(Tc2)>5
                tmp = fit(Tc2.takingoff+Tc2.landing, Tc2.ebird_q75_01, fittype('a*x'));
                fit_chg_ebird_y(i_c,i_y,1) = tmp.a;
                fit_chg_ebird_y(i_c,i_y,2) = sum(Tc2.ebird_n);
            end
            % figure; hold on; plot(Tc2.ts, Tc2.takingoff+Tc2.landing,'o');  plot(Tc2.ts, Tc2.ebird_q75_01,'-')
        end
    end
end


% All together
Tkt = Tkfit;
Tktoff = Tkt;
Tktoff.ts = Tktoff.ts+1;
Tc = innerjoin(Tkt,Tktoff,keys=["ts","name"],RightVariables=["ebird_q75","ebird_n"]);
Tc = renamevars(Tc,["ebird_q75_Tkt","ebird_q75_Tktoff"],["ebird_q75_0","ebird_q75_1"]);
Tc.ebird_q75_01 = Tc.ebird_q75_1-Tc.ebird_q75_0;
Tc.ebird_n = min([Tc.ebird_n_Tkt Tc.ebird_n_Tktoff],[],2);
% Tc=Tc(~isnan(Tc.ebird_n) & ~isnan(Tc.ebird_q75_01),:);

fit_chg_ebird_all=fit(Tc.takingoff+Tc.landing, Tc.ebird_q75_01, fittype('a*x'));
1./confint(fit_chg_ebird_all,.90);

%%

Tkt = Tk;
Tktoff = Tkt;
Tktoff.ts = Tktoff.ts+1;
Tc = innerjoin(Tkt,Tktoff,keys=["ts","name"],RightVariables=["ebird_q75","ebird_n"]);
Tc = renamevars(Tc,["ebird_q75_Tkt","ebird_q75_Tktoff"],["ebird_q75_0","ebird_q75_1"]);
Tc.ebird_q75_01 = Tc.ebird_q75_1-Tc.ebird_q75_0;
Tc.ebird_n = min([Tc.ebird_n_Tkt Tc.ebird_n_Tktoff],[],2);

figure; tiledlayout('flow','TileSpacing','none','Padding','tight')
nexttile; hold on ;
plot(Tc.ebird_n, abs(Tc.ebird_q75_01), '.k')
edges = [0:5:100 110:10:200 225:25:400 450:50:800];
Y = discretize(Tc.ebird_n,edges);
A=splitapply(@mean,abs(Tc.ebird_q75_01),Y);
plot(edges(1:end-1)+diff(edges),A,'.-r','MarkerSize',10)
grid on; box on; axis tight; ylim([0 40]); xlim([0 400])
ylabel("Daily difference in the count")
nexttile; hold on ;
ecdf(Tc.ebird_n); xlim([0 400])
grid on; box on;
xlabel("Number of checklists per day")
ylabel("Empirical cumulative distribution")


%% Figure

figure; box on; grid on; hold on;
scatter(Tc.takingoff+Tc.landing, Tc.ebird_q75_01,Tc.ebird_n/2,'ok','filled','MarkerFaceAlpha',.2)
plot(fit_chg_ebird_all); ylim([-30 30]);
% tmp=xlim; tmp=linspace(tmp(1), tmp(2),100);
% plot(tmp,1./predint(fit_chg_ebird_all,tmp,.5,'functional'))
ylabel("eBird change bird/checklist")
xlabel(" Weather Radar change birds/km^2")

if false
    y_lim=8;
    for i_c=1:height(cities)
        figure('position',[0 0 850 850]);  tiledlayout(3,1,'TileSpacing','tight','Padding','tight')

        Tkt=Tkfit;
        Tkt = Tkt(Tkt.ebird_n>thr_n&Tkt.name==cities.name(i_c),:);
        Tktoff=Tkt;
        Tktoff.ts = Tktoff.ts+1;
        Tc = innerjoin(Tkt,Tktoff,keys=["ts","name"],RightVariables="ebird_q75");
        Tc = renamevars(Tc,["ebird_q75_Tkt","ebird_q75_Tktoff"],["ebird_q75_0","ebird_q75_1"]);
        Tc.ebird_q75_01 = Tc.ebird_q75_1-Tc.ebird_q75_0;

        nexttile([2 1]); hold on; box on; grid on;
        x=(Tc.takingoff+Tc.landing);
        y=Tc.ebird_q75_01;
        x_ax = floor(min(x)):1:ceil(max(x));
        scatter(x,y,Tc.ebird_n*2,'.k',"AlphaData",.4)
        p11 = predint(fit_chg_ebird{i_c},x_ax,0.99,'functional','off');
        p = fill([x_ax fliplr(x_ax)]',[p11(:,1);flipud(p11(:,2))],'red', FaceAlpha = .3, EdgeColor = 'none');
        plot(fit_chg_ebird{i_c})
        xlim([x_ax(1) x_ax(end)])
        xlim([-40 40])
        ylim([-y_lim y_lim])
        xlabel('change on the ground = landing-takeoff'); ylabel('eBird q75'); title("Corr="+num2str(corr(x,y))); legend off

        nexttile; hold on; box on; grid on;
        %x=(Tc.takingoff-Tc.landing)/2;
        x=Tc.mvt;
        y=abs(Tc.ebird_q75_01);
        f=fit(x,y,"poly1",Weights=Tc.ebird_n);
        dx=.01;
        x_ax = (floor(min(x)/dx)*dx):dx:(ceil(max(x)/dx)*dx);
        p11 = predint(f,x_ax,0.99,'functional','off');
        scatter(x,y,Tc.ebird_n,'.k',"AlphaData",.4)
        p = fill([x_ax fliplr(x_ax)]',[p11(:,1);flipud(p11(:,2))],'red', FaceAlpha = .3, EdgeColor = 'none');
        plot(f)
        xlim([x_ax(1) x_ax(end)])
        xlabel('max(take-off,landing)'); ylabel('eBird q75'); title("Corr="+num2str(corr(x,y))); legend off
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
xlim([0 .5])
ylim([0 height(cities)+1])
yticks(1:height(cities))
yticklabels(cities.name)
xlabel('Correlation')

%
figure;
borders("states",'k');
scatter(cities.longitude, cities.latitude, c(:,3)/300, c(:,2), 'filled' ,"MarkerEdgeColor",[71 30 15]/255)
clim([0 .5]); colorbar; axis equal
colormap(crameri("batlow"))
axis([-110.9129  -62.9576   17.9394   59.6750])


t = table(...
    reshape(repmat(cities.name,1,numel(y)),[],1),...
    reshape(repelem(y,1,height(cities)),[],1),...
    reshape(fit_chg_ebird_y(:,:,1),[],1),...
    reshape(fit_chg_ebird_y(:,:,2),[],1),...
    VariableNames=["city", "year","a", "n"]  );
t.city = reordercats(categorical(t.city),cities.name);

figure('position',[0 0 800 350]); hold on;
% scatter(Tc.takingoff+Tc.landing./Tc.ebird_q75_01, categorical(Tc.name,cities.name),Tc.ebird_n/10,'ok','filled','MarkerFaceAlpha',.2)
scatter(1./t.a, t.city, t.n/100, 'k', 'filled','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
plot(fit_chg_ebird(:,2:3)', [categorical(cities.name,cities.name) categorical(cities.name,cities.name)]','-r')
scatter(fit_chg_ebird(:,1),categorical(cities.name,cities.name),c(:,3)/200,'r', 'filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','k')
grid on; box on;
xlabel("birds/km^2 equivalent to 1 bird/checklist")
xlim([-20 80])
%xlim([0 40])
xlim([-10 70])
xline(1./fit_chg_ebird_all.a,'r','LineWidth',2)
xline(quantile(repelem(1./t.a(t.n>0),t.n(t.n>0)),[.25 .75]),'--r')
set(gca,"ydir","reverse")



%% Coeff from seasonal match
param_match=nan(height(cities),numel(y),4);
for i_c=1:height(cities)
    for i_y=1:numel(y)
        Tkt = Tk(year(Tk.ts)==y(i_y)&Tk.name==cities.name(i_c),:);
        tmp = cumsum(-Tkt.takingoff-Tkt.landing,'omitnan');

        is_nan = ~isnan(Tkt.ebird_q75) & ~isnan(Tkt.ebird_n) & ~isnan(tmp);

        f = @(x) nansum(abs(Tkt.ebird_n(is_nan) .* (Tkt.ebird_q75(is_nan) - x(2).*tmp(is_nan) + x(1))))/sum(Tkt.ebird_n(is_nan));
        [param_match(i_c,i_y,1:2),fval]=fminsearch(f,[-50 0.05]);
        param_match(i_c,i_y,3)=fval;
        param_match(i_c,i_y,4)= sum(Tkt.ebird_n(is_nan));
    end
end

t = table(...
    reshape(repmat(cities.name,1,numel(y)),[],1),...
    reshape(repelem(y,1,height(cities)),[],1),...
    reshape(param_match(:,:,1),[],1),...
    reshape(param_match(:,:,2),[],1),...
    reshape(param_match(:,:,3),[],1),...
    reshape(param_match(:,:,4),[],1),...
    VariableNames=["city", "year","R", "slope", "fval", "n"]  );
t.city = reordercats(categorical(t.city),cities.name);
t = t(~isnan(t.fval),:);

figure('position',[0 0 800 350]); hold on;
scatter(1./t.slope, t.city, t.n/100, 'k', 'filled','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
grid on; box on; xlim([-10 70])
xlabel("birds/km^2 equivalent to 1 bird/checklist")
set(gca,"ydir","reverse")


%% Figure Appendix
ratio_city = min(max(0,fit_chg_ebird(:,1)),22);
ratio_city(:) = 16;

for i_c=1:height(cities)
    figure('position',[0 0 1100 850]); tiledlayout(5,3,'TileSpacing','tight','Padding','none');
    Tkt = Tk(Tk.name==cities.name(i_c),:);
    y_eb = nanstd(Tkt.ebird_q75)*[-2 2]+nansum(Tkt.ebird_q75.*Tkt.ebird_n) ./ nansum(Tkt.ebird_n);
    y_eb = round(y_eb/2)*2;
    y_wr = [0 diff(y_eb)*ratio_city(i_c)];
    for i_y=fliplr(y)
        Tkt = Tk(year(Tk.ts)==i_y&Tk.name==cities.name(i_c),:);

        if i_y==2021
            nexttile([2,2]);
        else
            nexttile;
        end
        hold on; box on; grid on;
        % patch([Tkt.ts+.5 ; flipud(Tkt.ts+.5)], [Tkt.ebird_q25 ; flipud(Tkt.ebird_q75)],cleBird2,"EdgeAlpha",0,"FaceAlpha",1)

        % plot(Tkt.ts+.5, Tkt.ebird_q25,'--',"color",cleBird)
        % plot(Tkt.ts+.5, Tkt.ebird_q75,'--',"color",cleBird)
        % a = cumsum(Tkt.ebird_n,'omitnan');
        % b = zeros(1,a(end));
        % b(a - Tkt.ebird_n +1)=1;
        % scatter(Tkt.ts(cumsum(b))+.5, vertcat(Tkt.ebird_d{:}),20,'filled',"MarkerFaceColor",cleBird,"MarkerFaceAlpha",.1)
        plot(Tkt.ts+.5, Tkt.ebird_q75,'-',"color",cleBird, LineWidth=2)
        scatter(Tkt.ts+.5, Tkt.ebird_q75, Tkt.ebird_n/2*(1+(i_y == 2021)),'filled',"MarkerFaceColor",cleBird)
        ylim(y_eb); y_ticks_eb = yticks;
        if i_y~=2021
            yticklabels('');
        end

        yyaxis right

        % bar(Tkt.ts,-Tkt.Fd_takingoff(:,i_c)-Fd_landing(:,i_c),1,'FaceColor','w');
        tmp = cumsum(-Tkt.takingoff-Tkt.landing,'omitnan');
        tmp = tmp-min(tmp)-range(tmp)/2 + y_wr(2)/2;
        bar(Tkt.ts,Tkt.takingoff,1,'FaceColor',clmap(end,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
        bar(Tkt.ts,-Tkt.landing,1,'FaceColor',clmap(1,:),FaceAlpha=1,EdgeColor=[.5 .5 .5])
        bar(Tkt.ts,min(Tkt.takingoff,-Tkt.landing),1,FaceAlpha=1,FaceColor=[.85 .85 .85],EdgeColor=[.5 .5 .5])
        plot(Tkt.ts+.5,tmp,'--k','LineWidth',3)

        ylim(y_wr);

        y_ticks_wr = interp1(y_eb,y_wr,y_ticks_eb);
        yticks(y_ticks_wr)
        if i_y~=2021
            yticklabels('');
        end

        xlim([ datetime(i_y,4,1) datetime(i_y,6,1)]);
        xticks(datetime(i_y,4:6,1)); xticklabels(string(datetime(i_y,4:6,1),"dd-MMM"));
        if i_y~=2021
            xticklabels([])
            txt = " "+num2str(i_y);
        else
            txt = " "+cities.name(i_c)+" | "+num2str(i_y);
        end
        text(datetime(i_y,4,1), y_wr(2),txt,HorizontalAlignment="left",VerticalAlignment="top",FontSize=24,FontWeight="bold",Color=[.2 .2 .2])
        set(gca,'XColor','k'); set(gca,'yColor','k');
    end
    % exportgraphics(gcf, "cumulative_"+num2str(i_c)+"_"+cities.name(i_c)+".png")
end

%% Abundance turnover

hexString = {'#e5ffca','#d0efb1','#c0e0ae','#b0d0ab','#a0c1a8','#8fb1a5','#7fa2a2','#6e929f','#5e829c','#4d7298'};
colt=hex2rgb(hexString);

ratio_city = fit_chg_ebird(:,1);
%ratio_city = min(max(7,ratio_city),22);
i_cc = find(ratio_city>0 & ratio_city<30);
% ratio_city(:)=16;

col=[interp1(linspace(0,1,numel(hexString))',colt(:,1),linspace(0,1,numel(i_cc)));...
    interp1(linspace(0,1,numel(hexString))',colt(:,2),linspace(0,1,numel(i_cc)));...
    interp1(linspace(0,1,numel(hexString))',colt(:,3),linspace(0,1,numel(i_cc)))]';

% ratio_city = sum(param_match(i_c,:,2) .* param_match(i_c,:,4)) ./ sum(param_match(i_c,:,4));

intervals = 0:.01:1;
histw = nan(numel(intervals),height(cities));
for i_c=1:height(cities)
    Tkt = Tk(Tk.name==cities.name(i_c),:);
    tmp1 = (abs(Tkt.takingoff)+abs(Tkt.landing))/2;
    tmp = tmp1/ratio_city(i_c)./Tkt.ebird_q75;
    histw(:,i_c) = histwc(tmp, tmp1, intervals);
end

figure('position',[0 0 850 200]); tiledlayout(3,4,'TileSpacing','tight','Padding','tight');
for i_c2=1:numel(i_cc)%height(cities)
    i_c=i_cc(i_c2);
    nexttile; hold on; box on;
    title(cities.name(i_c) + " | "+ num2str(round(ratio_city(i_c))))
    bar(intervals, histw(:,i_c), 1, "facecolor",col(i_c2,:))
    xlim([0 .5])
    grid on; xticklabels('');  yticks([]);
end

figure; b=bar(intervals, fliplr(histw(:,i_cc)),1,"stacked");
for k = 1:numel(b)
    b(k).FaceColor = col(size(col,1)-k+1,:);
end

intervals*histw(:,i_cc)./sum(histw(:,i_cc))

n = sum(histw(:,i_cc),'all');
mu = intervals*sum(histw(:,i_cc),2) ./ n;

sd = sqrt(((intervals-mu).^2*sum(histw(:,i_cc),2))/(n-1));



figure('position',[0 0 850 900]); tiledlayout(3,4,'TileSpacing','none','Padding','none');
for i_c=1:numel(i_cc)%height(cities)
    i_c=i_cc(i_c);
    nexttile; hold on; box on;


    Tkt = Tk(Tk.name==cities.name(i_c),:);
    Tkt.doy = day(Tkt.ts, 'doy');

    scatter(Tkt.doy, Tkt.ebird_q75, Tkt.ebird_n/4,'o','filled',"MarkerFaceColor",cleBird,'MarkerFaceAlpha',.2)
    for i_y=1:numel(y)
        %tmp = Tkt(year(Tkt.ts)== y(i_y),:);
        %plot(tmp.doy,smooth(tmp.ebird_q75),Color=max(0,1-(sum(tmp.ebird_n)/10000))*[1 1 1])
    end

    tmp = groupsummary(Tkt,"doy",@(x,y) sum(x.*y)/sum(y),{"ebird_q75", "ebird_n"});
    plot(tmp.doy,smooth(tmp.fun1_ebird_q75_ebird_n),'-',Color=cleBird, LineWidth=2)

    tmp = groupsummary(Tkt,"doy",@(x,y) quantile((abs(x)+abs(y))/2,[.05 .5 .95]),{"takingoff", "landing"});
    tmp2 = tmp.fun1_takingoff_landing./ratio_city(i_c);

    % bar(tmp.doy,tmp2(:,2),1,FaceColor=[.85 .85 .85],EdgeColor=[.5 .5 .5])
    errorbar(tmp.doy,tmp2(:,2), abs(tmp2(:,2)-tmp2(:,1)),abs(tmp2(:,2)-tmp2(:,3)),".",'Color',[.5 .5 .5],"LineStyle","none", "CapSize",0)

    axis tight;
    ylim([0 80]); grid on
    txt = cities.name(i_c) + " | "+ num2str(round(ratio_city(i_c)));
    text(91,80,txt,HorizontalAlignment="left",VerticalAlignment="top",FontSize=18,FontWeight="bold",Color=[.2 .2 .2])
    %title(txt)

    % figure; histogram((abs(Tkt.takingoff)+abs(Tkt.landing))/2/ratio_city(i_c)./Tkt.ebird_q75)

    xticks(day(datetime(2001,4:6,1),'doy'));
    xticklabels(["Apr" "May" "June"])
    % yticks(0:20:80)

end

%% Assess uncertainty
intervals = 0:.01:1;

n = 100;
histw2 = nan(numel(intervals),n);
histw = nan(numel(intervals),height(cities));
for i=1:n
    for i_c=1:height(cities)
        
        Tkt = Tk(Tk.name==cities.name(i_c),:);

        id = t.city==cities.name(i_c) & t.n>0;
        ratio = datasample(repelem(1./t.a(id),t.n(id)), 1);
        %[f,xi] = ksdensity(repelem(1./t.a(id),t.n(id)),"Bandwidth",10);
        %ratio = datasample(xi, 1, 'Weights', f);
        

        tmp1 = (abs(Tkt.takingoff));%+abs(Tkt.landing))/2;
        tmp = tmp1./ratio./Tkt.ebird_q75;
        histw(:,i_c) = histwc(tmp, tmp1, intervals);
    end
    histw2(:,i) = sum(histw(:,i_cc),2);
end
figure; hold on;
plot(intervals,smooth(quantile(histw2,.1,2)),'k')
plot(intervals,smooth(quantile(histw2,.5,2)),'k')
plot(intervals,smooth(quantile(histw2,.9,2)),'k')

ratio_city = fit_chg_ebird(:,1);
%ratio_city = min(max(7,ratio_city),22);
i_cc = find(ratio_city>0 & ratio_city<30);
figure; hold on;
ratio_list = [0, 16, 10, 20];

for i=1:numel(ratio_list)
    if ratio_list(i)==0
        ratio_city = fit_chg_ebird(:,1);
    else
        ratio_city(:) = ratio_list(i);
    end
    for i_c=1:height(cities)
        Tkt = Tk(Tk.name==cities.name(i_c),:);
        tmp1 = (abs(Tkt.takingoff));%+abs(Tkt.landing))/2;
        tmp = tmp1/ratio_city(i_c)./Tkt.ebird_q75;
        histw(:,i_c) = histwc(tmp, tmp1, intervals);
    end

    % plot(intervals, smooth(sum(histw(:,i_cc),2)), linewidth=2)
    plot(intervals, cumsum(sum(histw(:,i_cc),2))./sum(sum(histw(:,i_cc),2)), linewidth=2)
end
box on; grid on;
legend(string(ratio_list))



%% Seasonal curve year-year : Person's coefficient
Tkfit.doy = day(Tkfit.ts,'dayofyear');
% col=colormap('parula');
col = crameri("batlow");
clear c s10
plotit=false;
st=2010;
for i_c=1:height(cities)
    if plotit
        figure('position',[0 0 1600 900]); tiledlayout(2021-st+1,2021-st+1,'TileSpacing','none','Padding','none')
    end
    clear g
    for i_y1=st:2021
        Tkt1 = Tkfit(year(Tkfit.ts)==i_y1&Tkfit.name==cities.name(i_c)&Tkfit.ebird_n>0,:);
        s10{i_c}(i_y1-st+1) = sum(Tkt1.ebird_n<10);
        for i_y2=st:2021
            Tkt2 = Tkfit(year(Tkfit.ts)==i_y2&Tkfit.name==cities.name(i_c)&Tkfit.ebird_n>0,:);
            Tktm = outerjoin(Tkt1,Tkt2,LeftKeys=["doy","name"],RightKeys=["doy","name"],MergeKeys=true);

            tmp1 = normalize(cumsum(-Tktm.takingoff_Tkt1-Tktm.landing_Tkt1,'omitnan'),'zscore');
            tmp2 = normalize(Tktm.ebird_q75_Tkt2,'zscore');

            c{i_c}(i_y1-st+1,i_y2-st+1)=0;
            try
                c{i_c}(i_y1-st+1,i_y2-st+1)=corr(tmp1, tmp2,"rows",'complete');
            end

            %c(i,j)=sqrt(sum((zscore(tmp1(:,i))-zscore(tmp3(:,j))).^2));

            if plotit
                g(i_y1-st+1,i_y2-st+1)=nexttile; box on; hold on; % grid on;

                set(gca,'color',col(max(1,round(c{i_c}(i_y1-st+1,i_y2-st+1)*256)),:));
                if i_y1==i_y2
                    set(gca,'XColor','w'); set(gca,'yColor','w');
                else
                    set(gca,'XColor','k'); set(gca,'yColor','k');
                end
                xticklabels(''); yticklabels('');
                plot(tmp1,'--k');
                plot(tmp2,"color",cleBird);
            end
            % legend(num2str(c(i,j)))
        end
    end

    c{i_c}(isnan(c{i_c}))=0;
    M = matchpairs(1-abs(c{i_c}),1000);
    if plotit
        for i_m=1:height(M)
            set(g(M(i_m,1),M(i_m,2)),'XColor',"r"); set(g(M(i_m,1),M(i_m,2)),'yColor',"r");
            set(g(M(i_m,1),M(i_m,2)),'Linewidth',2);
        end
    end
    % exportgraphics(gcf, "figures/corr_"+num2str(i_c)+"_"+cities.name(i_c)+"_spring.png")
end



%% Seconal corr figure
% filter years with less than 5 days with less thans 10 checklists
id_valid = vertcat(s10{1:end-1})<5;
% filter city which have more than 5 valid year
id_city = find(sum(id_valid,2)>5);

figure('position',[0 0 800 350]); tiledlayout('flow','TileSpacing','tight','Padding','tight');
hold on; box on; grid on;
clear tmp
colormap(crameri("batlow"))
for i=1:numel(id_city)
    i_c = id_city(i);
    cs = abs(c{i_c}(id_valid(i_c,:),id_valid(i_c,:)));

    M = matchpairs(cs,1000);
    cmin = mean(cs(sub2ind(size(cs),M(:,1), M(:,2))));
    M = matchpairs(1-cs,1000);
    cmax = mean(cs(sub2ind(size(cs),M(:,1), M(:,2))));
    patch([cmin cmax nan]',[i i nan]',[cmin cmax nan]',[cmin cmax nan]', 'edgecolor', 'interp');

    %plot(cmin,i,'.',Color=cl(20,:),LineWidth=1, MarkerSize=40)
    %plot(cmax,i,'.',Color=cl(end-20,:),LineWidth=1, MarkerSize=40)
    scatter(cmin,i,20,cmin,'filled')
    scatter(cmax,i,20,cmax,'filled')

    % plot(mean(cs,'all'),i,'.',Color=cl(end/2+20,:),LineWidth=1, MarkerSize=40)

    plot(mean(diag(cs)),i,'.k', MarkerSize=15)
    tmp(i)=-(mean(diag(cs))-mean(cs(sub2ind(size(cs),M(:,1), M(:,2)))))/mean(cs(sub2ind(size(cs),M(:,1), M(:,2))));
end
set(gca,"ydir","reverse")
xlim([0 1])
ylim([0 numel(id_city)+1])
yticks(1:numel(id_city))
yticklabels(cities.name(id_city))
legend('Worse combinaison','Optimal combinaison','Random combinaison','True combinaison',FontSize=16,location="northoutside",Orientation="horizontal");
xlabel('Correlation')
% exportgraphics(gcf, "figures/corr_all_spring.png")
% exportgraphics(gcf, "figures_paper/figure_4c.eps")









%% Radar passage vs landing/takeoff

load('../BMM-US/data/density/inference-trans.mat')
vidT = vidTS;
vid = trans.f_inv(vidT);

Fd_dens = nan(numel(ts),height(cities));
idt = ismember(g.day,ts);

for i_c = 1:height(cities)
    i_c
    load("data/eBird/"+cities.name(i_c)+".mat")

    % i_r = radar.lat>=q_cond.latitude(1)-.2 & radar.lat<=q_cond.latitude(2)+.2 & radar.lon>=q_cond.longitude(1)-.2 & radar.lon<=q_cond.longitude(2)+.2;
    [distr,i_r]=mink(sqrt((radar.lat-mean(q_cond.latitude)).^2 + (radar.lon-mean(q_cond.longitude)).^2 )*111,5);
    w=1./distr'.* ~isnan(vid(:,i_r));
    w=w./sum(w,2);
    dens = sum(w.*vid(:,i_r),2,'omitnan');
    dens(all(isnan(vid(:,i_r)),2))=nan;

    tmp = splitapply(@nanmean,dens,g.day_id);
    Fd_dens(idt(idt>0),i_c) = tmp(idt>0);

    tmp = splitapply(@sum,sum(~isnan(vid(:,i_r)),2)>3,g.day_id)<30;
end

i_c=3;
clmap = crameri('berlin');
figure('position',[0 0 1650 1200]); tiledlayout('flow','TileSpacing','tight','Padding','tight'); set(gcf, 'color', 'k');
for i_y=2010:2021
    nexttile; box on; hold on
    bar(ts,Fd_takingoff(:,i_c),1,'FaceColor',clmap(end,:),'FaceAlpha',.5)
    bar(ts,-Fd_landing(:,i_c),1,'FaceColor',clmap(1,:),'FaceAlpha',.5)
    set(gca, 'color', 'k');
    set(gca,'XColor','w');
    set(gca,'yColor','w');
    yyaxis right;
    scatter(ts, Fd_dens(:,i_c), [],'o','filled');
    % scatter(g.day,tmp)
    set(gca,'yColor','w');
    xlim([datetime(i_y,4,1) datetime(i_y,6,1)])
end














%%

ts = datetime(2010,1,1):datetime(2022,1,1);
Fd_dens = nan(numel(ts),height(radar));
Fd_takingoff = nan(numel(ts),height(radar));
Fd_landing = nan(numel(ts),height(radar));
Fd_mvt = nan(numel(ts),height(radar));

radius=80/111;

for i_y=2010:2021
    i_y
    load(['../BMM-US/data/flow/est_' num2str(i_y) '.mat'])

    idt = find(ts==datetime(i_y,1,1))+(0:size(Fd.takingoff,3)-1);

    for i_r=1:height(radar)
        k = sqrt((g.LON-radar.lon(i_r)).^2 + (g.LAT-radar.lat(i_r)).^2)<=radius;
        k(isnan(Fd.takingoff(:,:,1))) = 0;
        k = k ./ sum(k,"all","omitnan");

        Fd_takingoff(idt,i_r)=squeeze(mean(Fd.takingoff./double(g.area) .* k,[1 2],'omitnan')); % bird/km^2
        Fd_landing(idt,i_r)=squeeze(mean(Fd.landing./double(g.area) .* k,[1 2],'omitnan')); % bird/km^2
        Fd_mvt(idt,i_r) =  squeeze(mean(MVT_day.*k,[1 2],'omitnan')); % bird/km (summed over the night - averaged over spaced
    end
end


%%
id= Fd_mvt(:)>nanmean(Fd_mvt(:));

corr(Fd_mvt(id),Fd_takingoff(id))


corr(reshape(Fd_mvt(idt{1},:),[],1), reshape(Fd_takingoff(idt{1},:),[],1))
corr(Fd_mvt(idt{1},:),Fd_takingoff(idt{1},:)+Fd_landing(idt{1},:))

figure; plot(Fd_mvt(id),abs(Fd_takingoff(id)+Fd_landing(id)),'.k')


%% spatial
clear idt
idt{1} = month(ts)>=4 & month(ts)<=5;
idt{2} = month(ts)>=9 & month(ts)<=10;
season=["Spring" "Autumn"];

figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')
for i=1:2
    nexttile;
    scatter(radar.lon, radar.lat, 200, nanmean(Fd_takingoff(idt{i},:),1),'filled');
    colorbar; axis equal tight; a_xis=axis();  borders("states",'k'); axis(a_xis); title("Takingoff (bird/km^2/night)")
    nexttile;
    scatter(radar.lon, radar.lat, 200, nanmean(-Fd_landing(idt{i},:),1),'filled');
    colorbar; axis equal tight;a_xis=axis();  borders("states",'k'); axis(a_xis); title("Landing (bird/km^2/night)")
    nexttile;
    scatter(radar.lon, radar.lat, 200, nanmean(Fd_takingoff(idt{i},:)+Fd_landing(idt{i},:),1),'filled');
    colorbar; axis equal tight;a_xis=axis();  borders("states",'k'); axis(a_xis); title("Change on the ground (bird/km^2/night)")
    nexttile;
    scatter(radar.lon, radar.lat, 200, nanmean(Fd_mvt(idt{i},:),1),'filled');
    colorbar; axis equal tight;a_xis=axis();  borders("states",'k'); axis(a_xis); title("Movement (bird/km/night)")
end

figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')
for i=1:2
    nexttile;
    tmp = nanmean((Fd_takingoff(idt{i},:)-Fd_landing(idt{i},:))/2,1)./nanmean(Fd_mvt(idt{i},:),1);
    scatter(radar.lon, radar.lat, 200, tmp,'filled');
    colorbar; axis equal tight;a_xis=axis();  borders("states",'k'); axis(a_xis); title("Stopover (bird/km^2/night)  / Passage (bird/km/night) "+ season(i))
    clim([0.1 .4])
end



tmp_1 = min(cat(3,Fd_takingoff,-Fd_landing),[],3);
tmp_2 = abs(Fd_takingoff+Fd_landing);

X=[tmp_1(:),tmp_2(:)];
X(X>500)=nan;
ksdensity(X)

figure; hold on;
plot(Fd_mvt,Fd_takingoff,'.r')
plot(Fd_mvt,Fd_landing,'.b')

figure; hold on;
plot(Fd_takingoff,Fd_landing,'.k')

figure; hold on;
plot(Fd_mvt,Fd_takingoff+Fd_landing,'.k')



%%

cgt = nan(numel(ts),2);
cgl = nan(numel(ts),2);
w = nan(numel(ts),1);

for i_y=2010:2021
    i_y
    load(['../BMM-US/data/flow/est_' num2str(i_y) '.mat'])

    Fd.takingoff(isnan(Fd.takingoff))=0;
    Fd.landing(isnan(Fd.landing))=0;

    idt = find(ts==datetime(i_y,1,1))+(0:size(Fd.takingoff,3)-1);

    for i_v=1:numel(idt)
        props = regionprops(true(size(Fd.takingoff(:,:,1))), Fd.takingoff(:,:,i_v), 'WeightedCentroid');
        id = floor(props.WeightedCentroid);
        if ~any(isnan(id))
            cgt(idt(i_v),:) = [g.lon(id(1)) g.lat(id(2))] + (props.WeightedCentroid-id).*diff(g.lat(1:2));
        end
        props = regionprops(true(size(Fd.takingoff(:,:,1))), Fd.landing(:,:,i_v), 'WeightedCentroid');
        id = floor(props.WeightedCentroid);
        if ~any(isnan(id))
            cgl(idt(i_v),:) = [g.lon(id(1)) g.lat(id(2))] + (props.WeightedCentroid-id).*diff(g.lat(1:2));
        end
    end
    w(idt)=squeeze(nansum(MVT_day,[1 2]));
end



figure; tiledlayout('flow','TileSpacing','tight','Padding','tight')%
% scatter(cgt(:,1),cgt(:,2),w/3e6,1:(size(tmp,3)-1))
sm_p =1;
for i_y=2010:2021
    nexttile; hold on; borders('continental us','k')
    id = year(ts)==i_y & (month(ts)==4 | month(ts)==5);
    plot(smooth(cgt(id,1),sm_p), smooth(cgt(id,2),sm_p),'-r')
    scatter(smooth(cgt(id,1),sm_p),smooth(cgt(id,2),sm_p),w(id)/100000,datenum(ts(id)),'filled')

    plot(smooth(cgl(id,1),sm_p), smooth(cgl(id,2),sm_p),'-b')
    scatter(smooth(cgl(id,1),sm_p),smooth(cgl(id,2),sm_p),w(id)/100000,datenum(ts(id)),'filled')
    axis equal; axis([-104 -92 35 42]);
end



distgc = deg2km(distance(cgt(:,2),cgt(:,1),cgl(:,2),cgl(:,1)));


figure;
id = ~isnan(w) & ~isnan(distgc);
scatter(day(ts(id),'dayofyear'),sign(cgl(id,2)-cgt(id,2)).*distgc(id),w(id)/100000,year(ts(id)),'filled')
ylabel("Distance (S->N)"); box on; grid on;
yline(0,'--k',LineWidth=2)
xlim([50 340]); ylim([-300 300])


%%
tmp_1 = squeeze(sum(min(cat(4,Fd.takingoff,-Fd.landing),[],4),[1 2],"omitnan"));

tmp_2 = squeeze(sum(abs(Fd.takingoff+Fd.landing),[1 2],"omitnan"));

figure; hold on;
bar(tmp_1, facealpha=.5)
bar(tmp_2, facealpha=.5);
legend("min(takeoff,landing)", "abs(take-off-landing)");
grid on; box on;



Fd.takingoff(:,:,100),  Fd.takingoff(:,:,100)






