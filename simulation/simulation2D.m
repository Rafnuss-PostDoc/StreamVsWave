% Set grid size
nx=120;
ny=40;
%nz=10;
nt=12*60/5;

% Generate Departure
addpath('/Users/raphael/Documents/Github/FastGaussianSimulation')
sim.s=[ny nx nt];
sim.n = 1;
covar(1).model = 'k-bessel';
covar(1).range = [4 10 5];
covar(1).azimuth = [0 0];
covar(1).var = .8;
covar(1).alpha = 5;
covar(2).model = 'k-bessel';
covar(2).range = [1 1 1];
covar(2).azimuth = [0 0];
covar(2).var = 0.01;
covar(2).alpha = 5;
D = FGS(sim,covar);
D = D{1};
D(D<0) = 0;
wt = 1:nt;
x=(wt-mean(wt))/3;
wt_f = (2*(1-exp(x)./(1+exp(x))))-1;
% plot(wt,wt_f)
D = D .* reshape(wt_f,1,1,[]);
D(D<0) = 0;

% Prob of landing 
wt_f2 = exp((wt-wt(end))/2);
L = ones(size(D)) .* reshape(wt_f2,1,1,[]);

% Test landing departure
figure('position',[0 0 1600 600]);
for i_t=1:2:(nt/2)
    imagesc(D(:,:,i_t)); title(i_t/nt)
    colorbar; axis tight equal; caxis([0 2])
    drawnow;pause(.05);
end

% Generate speed
U = -50*ones(ny,nx,nt);%+ normrnd(0,.01,ny,nx,nt);
V = 10*ones(ny,nx,nt);%+ normrnd(0,.01,ny,nx,nt);
% W = 1*ones(nx,ny,nz,nt).*reshape(wt_f,1,1,1,[]);

% dd=4;
% [Y,X] = meshgrid(1:nx,1:ny);
% figure('position',[0 0 1600 600]);tiledlayout(3,1,'TileSpacing','none','Padding','none')
% nexttile;quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),U(1:dd:end,1:dd:end,1),V(1:dd:end,1:dd:end,1))
% nexttile;quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),U(1:dd:end,1:dd:end,end/2),V(1:dd:end,1:dd:end,end/2))
% nexttile;quiver(X(1:dd:end,1:dd:end),Y(1:dd:end,1:dd:end),U(1:dd:end,1:dd:end,end),V(1:dd:end,1:dd:end,end))

%% Forward model
DS=nan(ny,nx,nt); % density
A=nan(ny,nx,nt); % arrival
% initialize to zero
DS(:,:,1)=0;
dx=25; dy=25; dt=5/60;

for i_t=1:nt-1
    rho = DS(:,:,i_t);
    vx = U(:,:,i_t);
    vy = V(:,:,i_t);
    
%     rho=[0:10 10:-1:0];
%     vx=ones(size(rho));
    
    % *Compute the flux* $\Phi =\mathit{\mathbf{v}}\rho$ *at +/- 1/2 grid cell.*
    % First, add a nan layer in lat or lon direction
    Phiy_pad = padarray(rho .* vy, [1 0 0], nan); % bird * km/h
    Phix_pad = padarray(rho .* vx, [0 1 0], nan);

    % and then, compute the flux at +/- 1/2 even if either previous or next cells
    % is nan.
    Phiy_h = movmean(Phiy_pad,[0 1],1,'omitnan','Endpoints','discard');
    Phix_h = movmean(Phix_pad,[0 1],2,'omitnan','Endpoints','discard');

    % *Compute the delta flux / delta distance* $\frac{\Delta \Phi }{\Delta \left(\textrm{lat},\textrm{lon}\right)}$
    % First, add 0 padding for the outer zone (allows to compute boundary cell)
    Phiy_h_0=padarray(Phiy_h,[1 1],0);
    Phix_h_0=padarray(Phix_h,[1 1],0);

    % Finally, compute the delta flux over delta distance.
    dPhidy = diff(Phiy_h_0,1,1)/dy; % bird * km/h /km -> bird/h
    dPhidx = diff(Phix_h_0,1,2)/dx;
    
    % *Compute the variation of bird in/out each cell*
    F = (dPhidy + dPhidx ).*dt; % bird/h * hr -> bird
    
    % F comprise both the intral change of bird and the bird going out of the area
    % Inner flux: remove the boundary cells
    Fin = F;
    Fin([1 end],:) = 0; Fin(:,[1 end]) = 0;

    % Outer Flux. Note that Philat_h_0 and Philat_h_0 are 0 for the cell outside
    % the domain, such at dPhilatdlat is equal to the flux at the outer cell (with
    % the sign corresponding to the direction).
    % Fout = F;
    % Fout(repmat(~gext.mask_water,1,1,gext.nt)) = 0;

    % Compute next time step
     rho1 = rho + Fin(2:end-1,2:end-1) + D(:,:,i_t);
     A(:,:,i_t) = rho1 .* L(:,:,i_t);
     DS(:,:,i_t+1) = rho1-A(:,:,i_t);
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

figure('position',[0 0 1600 600]);
for i_t=1:2:nt
    imagesc(DS(:,:,i_t)); title(i_t*dt)
    colorbar; axis tight equal; caxis([0 max(DS(:))])
    drawnow;pause(.01);
end

%% sampling density
DS_known = nan(size(DS));
DS_known(10:20,20:30,:) = DS(10:20,20:30,:);
DS_known(30:40,100:110,:) = DS(30:40,100:110,:);

figure; imagesc(DS_known(:,:,10));

%% What you have access to:
% - U,V
% - DS(:,:,1)=0; DS(:,:,end)=0;
% - DS_known
% - D+L is smooth

