% Set grid size
dx=20; % km
dt=1/60; % hr
nx=2000/dx; % km -> nb cell
nt=10/dt; % hr -> nb of cell

x=(1:nx)*dx;
t=(1:nt)*dt;

clmap = crameri('berlin');

[T, X] = meshgrid(1:nt, 1:nx);

%% Generate data
% Departure
r = 500;
fx = normpdf((1:nx)',.5*nx, r/dx);
fx = 0*max(fx)+fx;

% Landing
if false
    u=((1:nt)-nt/2)/(nt/2);
    u=(u+1).^.65/2.^.65*2-1;
    ft = (1-abs(u).^2);%normpdf((1:nt),.3*nt,nt/4);
    ft = ft ./sum(ft);
    d =  fx .* ft;
    d = 5000.*d./sum(d(:));
    l=1-d;
else
    ft = normpdf((1:nt),0,nt/10);
    ft = ft ./sum(ft);
    d =  fx .* ft;
    d = d./sum(d(:));

    ft2 = normpdf((1:nt),1*nt,nt/10);
    ft2 = ft2 ./sum(ft2);
    l =  ones(size(fx)) .* ft2;
    l = l./max(l(:));
end

% Generate speed
% Uniform
U = -34*ones(nx,nt);% km/hr
%U = U + fx*600; % faster at optimal condition
abs(U(1,1)*dt/dx)

assert(all(d(:)<=1))
assert(all(l(:)<=1))

if false
    figure('position',[0 0 600 200]); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    nexttile; plot(x,fx); xlabel('Space (km)'); yticks('')
    nexttile; plot(t,ft); xlabel('Time (hours)');yticks('')
    nexttile; imagesc(x,t,d); xlabel('Space'); ylabel('Time')
end

if false
    figure('position',[0 0 600 200]); tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    nexttile; plot(x,fx); xlabel('Space (km)'); yticks('')
    nexttile; plot(t,ft); xlabel('Time (hours)');yticks('')
    nexttile; plot(t,ft2); xlabel('Time (hours)');yticks('')
end

% barrier
% wt_f3 = min(1,exp(((1:nx)-nx*.9)/2));
% D = D .* (1-wt_f3');
% L = min(1,L+wt_f3');

% Test landing departure
if false
    figure('position',[0 0 1600 600]);
    for i_t=1:2:(nt/2)
        plot(d(:,i_t)); title(i_t/nt)
        axis tight; ylim([0 max(d(:))])
        drawnow;pause(.05);
        xlabel('Space'); ylabel('Density')
    end
end

%% Forward model
A=nan(nx,nt); % air
L=nan(nx,nt); % landing
D=nan(nx,nt); % departure
G = nan(nx,nt); % ground

% Initialize
A(:,1) = 0; % no bird in the air
G(:,1) = 1; % 1 bird on the gorund bird/km

for i_t=1:nt-1

    
    % *Compute the flux* $\Phi =\mathit{\mathbf{v}}\rho$ *at +/- 1/2 grid cell.*
    % First, add a nan layer in lat or lon direction
    %Phix_pad = padarray(rho .* vx, [1 0], 0);% bird * km/h
    Phix_pad = padarray(A(:,i_t) .* U(:,i_t), [1 0], nan);
    Phix_pad = fillmissing(Phix_pad,'pchip');
    Phix_pad(1) = max(Phix_pad(1),Phix_pad(2));

    % and then, compute the flux at +/- 1/2 even if either previous or next cells
    % is nan.
    Phix_h = movmean(Phix_pad,[0 1],'omitnan','Endpoints','discard');

    % *Compute the delta flux / delta distance* $\frac{\Delta \Phi }{\Delta \left(\textrm{lat},\textrm{lon}\right)}$
    % First, add 0 padding for the outer zone (allows to compute boundary cell)
    Phix_h_0=padarray(Phix_h,1,0);

    % Finally, compute the delta flux over delta distance.
    dPhidx = diff(Phix_h_0)/dx;% bird * km/h /km -> bird/h
    
    % *Compute the variation of bird in/out each cell*
    F = dPhidx.*dt; % bird/h * hr -> bird
    
    % F comprise both the intral change of bird and the bird going out of the area
    % Inner flux: remove the boundary cells
    Fin = F;
    Fin([1 end]) = 0; 

    % Outer Flux. Note that Philat_h_0 and Philat_h_0 are 0 for the cell outside
    % the domain, such at dPhilatdlat is equal to the flux at the outer cell (with
    % the sign corresponding to the direction).
    % Fout = F;
    % Fout(repmat(~gext.mask_water,1,1,gext.nt)) = 0;

    % Compute next time step
    A(:,i_t+1) = A(:,i_t)+Fin(2:end-1);
    %assert(all(A(:,i_t+1)>=0))

    D(:,i_t) = G(:,i_t).*d(:,i_t);
    L(:,i_t) = A(:,i_t+1).*l(:,i_t);
    tmp = D(:,i_t)-L(:,i_t);
    D(:,i_t) = tmp; D(tmp<0,i_t)=0;
    L(:,i_t) = -tmp; L(tmp>0,i_t)=0;

    A(:,i_t+1) = A(:,i_t+1)+tmp;
    G(:,i_t+1) = G(:,i_t)-tmp;
    %assert(all(A(:,i_t+1)>=0))
    %assert(all(G(:,i_t+1)>=0))
    

     % illustration
     if false
         figure('position',[0 0 2400 1000]); tiledlayout(5,1,'TileSpacing','tight','Padding','tight');
         nexttile; hold on; plot(A(:,i_t));plot(A(:,i_t+1)); title('Air');  
         nexttile; hold on; plot(G(:,i_t));plot(G(:,i_t+1)); title('Ground');  
         nexttile; plot(D(:,i_t)); title('Departure');  
         nexttile; plot(L(:,i_t)); title('Landing'); 
         nexttile; plot(Fin(2:end-1)); title('Movment');
         keyboard
     end
end

 %% Vizalize

clord = colororder;
% figure; imagesc(DS); xlabel('time'); ylabel('space')
% 
% figure;  hold on;
% plot(-nansum(A,2)); 
% plot(sum(D,2)); 
% plot(sum(D,2)); 
% xlabel('space'); ylabel('time')
% 
% 
% figure;  hold on;
x_radar = round(nx/2+[-1 0 1]*200/dx);

% plot(DS(x_radar,:)'); 
% xlabel('time'); ylabel('density')

%

figure('position',[0 0 1000 550]);tiledlayout(3,1,'TileSpacing','tight','Padding','tight');
set(gcf, 'color', 'none'); 

ax2=nexttile; hold on; ax2.Color="k";ax2.XColor="w"; ax2.YColor="w"; xticks([]); yticks([])
p1 = plot(A(:,end),'w','LineWidth',2);  ylim([0 max(A(:))])
p4 = scatter(x_radar,A(x_radar,end)',100,clord(1:3,:),'filled');
ylim([0 max(A(:))]); xlim([0 nx]); ylabel('Flying','FontSize',16); xlabel('Space','FontSize',16); box on; 

ax1 = nexttile; hold on; ax1.Color="k";ax1.XColor="w"; ax1.YColor="w"; xticks([]); yticks([])
plot(G(:,1),'w')
p2 = plot(sum(L,2,'omitnan'),'color',clmap(1,:),'LineWidth',2); 
p3 = plot(sum(-D,2,'omitnan'),'color',clmap(end,:),'LineWidth',2);  
p7 = plot(G(:,i_t),'--w','LineWidth',2);  
%ylim([0 max(G(:,1)+sum(L,2,'omitnan'))])
ylim([min(G(:,1)-sum(D,2,'omitnan')) max(G(:,1)+sum(L,2,'omitnan'))])
xlabel('Space','FontSize',16); ylabel('Departing(-) | Landing(+)','FontSize',16); box on;

ax3=nexttile; hold on; ax3.Color="k";ax3.XColor="w"; ax3.YColor="w"; xticks([]); yticks([])
xlim([0 nt]); ylim([0 max(max(A(x_radar,:)))])
ylabel('Flying','FontSize',16); xlabel('Time','FontSize',16)
p5 = plot(A(x_radar,:)','LineWidth',2); 
p6 = scatter(ones(1,3),A(x_radar,1)',100,clord(1:3,:),'filled');box on;

% video
for i_t=1:round(nt/100):nt
    p1.YData = A(:,i_t);
    p4.YData = max(0,A(x_radar,i_t));

    p2.YData =  G(:,1)+sum(L(:,1:i_t),2); 
    p3.YData =  G(:,1)+sum(-D(:,1:i_t),2); 
    p7.YData =  G(:,i_t); 
    p5(1).YData = A(x_radar(1),1:i_t)';
    p5(2).YData = A(x_radar(2),1:i_t)';
    p5(3).YData = A(x_radar(3),1:i_t)';
    p6.XData = i_t*ones(1,3); 
    p6.YData = max(0,A(x_radar,i_t)); 
    title(ax1,"Time: "+num2str(round(i_t/nt*100))+" %");
    drawnow;pause(.01);
    %exportgraphics(gcf,'simulation/flow_1.gif','Append',true,'BackgroundColor','k');
end
%exportgraphics(gcf,'simulation/flow_1.png','BackgroundColor','k')


% delete(p6)
% delete(p1)
% delete(p4)
% axes(ax2)
% for i_t=round([.1 .5 .9]*nt)
%     plot(DS(:,i_t),'w','LineWidth',2);  ylim([0 max(DS(:))])
%     scatter(x_radar,DS(x_radar,i_t)',100,clord(1:3,:),'filled');
% end

%% 

figure('position',[0 0 1000 550]);tiledlayout(2,1,'TileSpacing','tight','Padding','tight');
nexttile;hold on; xticks([]); yticks([])
plot(G(:,1),'k')
p2 = plot(G(:,1)+sum(L,2,'omitnan'),'color',clmap(1,:),'LineWidth',2); 
p3 = plot(G(:,1)+sum(-D,2,'omitnan'),'color',clmap(end,:),'LineWidth',2);  
p7 = plot(G(:,i_t),'--k','LineWidth',2);  
%ylim([0 max(G(:,1)+sum(L,2,'omitnan'))])
ylim([min(G(:,1)-sum(D,2,'omitnan')) max(G(:,1)+sum(L,2,'omitnan'))])
xlabel('Space','FontSize',16); ylabel('Departing(-) | Landing(+)','FontSize',16); box on;
exportgraphics(gcf,'figures_paper/figure_1.eps')
