% Set grid size
nx=120;
nt=12*60/5*16;
clmap = crameri('berlin');

%% Generate Departure
addpath('/Users/raphael/Documents/Github/FastGaussianSimulation')
sim.s=[nx nt];
sim.n = 1;
covar={};
covar(1).model = 'k-bessel';
covar(1).range = [10 6];
covar(1).azimuth = [0];
covar(1).var = .8;
covar(1).alpha = 5;
covar(2).model = 'k-bessel';
covar(2).range = [1 1];
covar(2).azimuth = [0];
covar(2).var = 0.01;
covar(2).alpha = 5;
D = FGS(sim,covar);
D = D{1};
D(D<0) = 0;
% imagesc(D)
wt = 1:nt;
x=(wt-mean(wt))/3;
s=4;
wt_f = (2*(1-exp(x/s)./(1+exp(x/s))))-1;
% plot(wt,wt_f)
D = D .* wt_f;
D(D<0) = 0;

%% Gaussian 
[tmp1, tmp2]=meshgrid(1:nt, 1:nx);
% Wave
% D = 2000*reshape(mvnpdf([tmp1(:) tmp2(:)],[.1 .25].*[nt nx],[6 2].*[nt nx]),nx,nt);

%D = 2000*reshape(mvnpdf([tmp1(:) tmp2(:)],[20 30],[1600 400])+...
%    1*mvnpdf([tmp1(:) tmp2(:)],[60 50],[8000 10000]),nx,nt);

% Flow
D = 20000*reshape(mvnpdf([tmp1(:) tmp2(:)],[.3 .4].*[nt nx],[4000 10].*[nt nx]),nx,nt);

% D = D.*(1+randn(size(D))*.05);

% figure; imagesc(D); ylabel('Space'); xlabel('Time')

% Prob of landing 
wt = 1:nt;

% wave
% wt_f2 = exp((wt-nt)/nt*40);


% flow
%x=(wt- (nt*1.2))/3; s=48;
%wt_f2 = 1-(1-exp(x/s)./(1+exp(x/s)));
wt_f2 = 2*normpdf((1:nx)',.6*nx,5*nx);
wt_f2 = wt_f2+normpdf((1:nt),nt,1*nt) + exp(((1:nt)-1.4*nt)/nt*10);%
imagesc(wt_f2)

% figure; plot(wt_f2)
L = ones(size(D)) .* wt_f2;

% figure; imagesc(log(L)); ylabel('Space'); xlabel('Time');caxis([-6.5 -5])

% barrier
% wt_f3 = min(1,exp(((1:nx)-nx*.9)/2));
% D = D .* (1-wt_f3');
% L = min(1,L+wt_f3');

% Test landing departure
% figure('position',[0 0 1600 600]);
% for i_t=1:2:(nt/2)
%     plot(D(:,i_t)); title(i_t/nt)
%     axis tight; ylim([0 max(D(:))])
%     drawnow;pause(.05);
%     xlabel('Space'); ylabel('Density')
% end

% Generate speed
U = -800*nx/nt*ones(nx,nt);%+ normrnd(0,.01,ny,nx,nt);
% U = -400*nx/nt*ones(nx,nt);
% W = 1*ones(nx,ny,nz,nt).*reshape(wt_f,1,1,1,[]);

% dd=4;
% [Y,X] = meshgrid(1:nx,1:ny);
% figure('position',[0 0 1600 600]);tiledlayout(3,1,'TileSpacing','none','Padding','none')
% nexttile;quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),U(1:dd:end,1:dd:end,1),V(1:dd:end,1:dd:end,1))
% nexttile;quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),U(1:dd:end,1:dd:end,end/2),V(1:dd:end,1:dd:end,end/2))
% nexttile;quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),U(1:dd:end,1:dd:end,end),V(1:dd:end,1:dd:end,end))

%% Forward model
DS=nan(nx,nt); % density
A=nan(nx,nt); % arrival
% initialize to zero
DS(:,1)=0;
dx=25; dy=25; dt=5/60/4;

for i_t=1:nt-1
    rho = DS(:,i_t);
    vx = U(:,i_t);
    
%     rho=[0:10 10:-1:0];
%     vx=ones(size(rho));
    
    % *Compute the flux* $\Phi =\mathit{\mathbf{v}}\rho$ *at +/- 1/2 grid cell.*
    % First, add a nan layer in lat or lon direction
    %Phix_pad = padarray(rho .* vx, [1 0], 0);% bird * km/h
    Phix_pad = padarray(rho .* vx, [1 0], nan);
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
    F = (dPhidx).*dt; % bird/h * hr -> bird
    
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
     rho1 = rho + Fin(2:end-1) + D(:,i_t);
     A(:,i_t) = rho1 .* L(:,i_t);
     DS(:,i_t+1) = rho1-A(:,i_t);
     % illustration
%      figure('position',[0 0 2400 1000]); tiledlayout('flow','TileSpacing','none','Padding','none'); colorbar;
%      nexttile; imagesc(rho); title('density current timestep');  colorbar;
%      nexttile; imagesc(Fin(2:end-1,2:end-1)); title('internal change'); colorbar;
%      nexttile; imagesc(rho+Fin(2:end-1,2:end-1)); title('internal change'); colorbar;
%      nexttile; imagesc(D(:,:,i_t)); title('departure'); colorbar;
%      % nexttile; quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),vx(1:dd:end,1:dd:end),vy(1:dd:end,1:dd:end))
%      nexttile; imagesc(A(:,:,i_t)); title('Arrival'); colorbar;
%      nexttile; imagesc(DS(:,:,i_t+1)); title('density next timestep'); colorbar;
%      keyboard
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
x_radar = round(nx*[.25 .5 .75]);
% plot(DS(x_radar,:)'); 
% xlabel('time'); ylabel('density')

%

figure('position',[0 0 1600 900]);tiledlayout(3,1,'TileSpacing','tight','Padding','tight');
set(gcf, 'color', 'none'); 

ax1 = nexttile; box on;hold on; ax1.Color="k";ax1.XColor="w"; ax1.YColor="w"; xticks([]); yticks([])
yline(0,'w')
p2 = plot(nansum(A(:,1),2),'color',clmap(1,:),'LineWidth',2); 
p3 = plot(-sum(D(:,1),2),'color',clmap(end,:),'LineWidth',2);  
p7 = plot(nansum(A(:,1),2)-sum(D(:,1),2),'--w','LineWidth',2);  
ylim([-max(sum(D,2)) max(nansum(A,2))])
xlabel('Space','FontSize',16); ylabel('Departing(-) | Landing(+)','FontSize',16)

ax2=nexttile; box on; hold on; ax2.Color="k";ax2.XColor="w"; ax2.YColor="w"; xticks([]); yticks([])
p1 = plot(DS(:,1),'w','LineWidth',2);  ylim([0 max(DS(:))])
p4 = scatter(x_radar,DS(x_radar,1)',100,clord(1:3,:),'filled');
ylim([0 max(DS(:))]); xlim([0 nx]); ylabel('Flying','FontSize',16); xlabel('Space','FontSize',16)

ax3=nexttile; box on; hold on; ax3.Color="k";ax3.XColor="w"; ax3.YColor="w"; xticks([]); yticks([])
xlim([0 nt]); ylim([0 max(max(DS(x_radar,:)))])
ylabel('Flying','FontSize',16); xlabel('Time','FontSize',16)
p5 = plot(DS(x_radar,:)','LineWidth',2); 
p6 = scatter(ones(1,3),DS(x_radar,1)',100,clord(1:3,:),'filled');

for i_t=1:round(nt/100):nt
    p1.YData = DS(:,i_t);
    p4.YData = max(0,DS(x_radar,i_t));
    p2.YData = nansum(A(:,1:i_t),2); 
    p3.YData = -sum(D(:,1:i_t),2); 
    p7.YData = nansum(A(:,1:i_t),2)-sum(D(:,1:i_t),2); 
    p5(1).YData = DS(x_radar(1),1:i_t)';
    p5(2).YData = DS(x_radar(2),1:i_t)';
    p5(3).YData = DS(x_radar(3),1:i_t)';
    p6.XData = i_t*ones(1,3); 
    p6.YData = max(0,DS(x_radar,i_t)); 
    title(ax1,"Time: "+num2str(round(i_t/nt*100))+" %");
    drawnow;pause(.01);
    exportgraphics(gcf,'simulation/flow_1.gif','Append',true,'BackgroundColor','k');
end
exportgraphics(gcf,'simulation/flow_1.png','BackgroundColor','k')



