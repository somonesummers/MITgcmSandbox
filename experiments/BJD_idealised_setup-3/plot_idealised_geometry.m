% plot_idealised_geometry

clear
close all
clc

%% Grid
% Cell resolution
deltaX = 500;
deltaY = 500;
deltaZ = 10;

% Dimensions of grid
nx=100; % 50 km long
ny=12; % 10 km wide (plus 500 m walls)
nz=50; % 500 m deep

% x scale
delx = zeros(1,nx); 
delx(:) = deltaX;

% y scale
dely = zeros(1,ny);
dely(:) = deltaY;

% Distances for plotting
Xdist = cumsum(delx./1000);
Ydist = cumsum(dely./1000);

% Sill
sill_depth = 100; % metres
h = sill_depth/(nz*deltaZ);

% Bathymetry
bathymetry = zeros(nx,ny); % pre-allocate
bathymetry(:) = -deltaZ*nz; % uniformly 500 m deep (negative)
bathymetry(1,:) = 0; % barrier behind virtual glacier front (at western end of domain)
bathymetry(:,[1 end]) = 0; % fjord walls

% Bathymetry with 100 sill
bathymetry_sill = zeros(nx,ny); % pre-allocate
bathymetry_sill(:) = -deltaZ*nz; % uniformly 500 m deep (negative)
%   - create gaussian sill
a1 = (deltaZ*nz)-((deltaZ*nz)*h); % amplitude - sill is 25% of fjord depth below the fjord surface
b1 = 10; % location of peak along sill (relates to x below)
c1 = 3; % relates to width of peak
x = 1:20; % number of cells wide
sill = a1.*exp(-((x-b1)./c1).^2);
bathymetry_sill(nx-30:nx-30+length(sill)-1,:) = bathymetry_sill(nx-30:nx-30+length(sill)-1,:) + repmat(sill(:),[1,size(bathymetry_sill,2)]);
% Fjord walls
bathymetry_sill(1,:) = 0; % barrier behind virtual glacier front (at western end of domain)
bathymetry_sill(:,[1 end]) = 0; % fjord walls

%% Initial and open boundary conditions
% Vertical cell spacing
zprof = -((0.5*deltaZ):deltaZ:((nz*deltaZ)-(0.5*deltaZ)));
delz = zeros(1,nz);
delz(:) = deltaZ;

% pre-allocate
temp = zeros([length(zprof),5]).*NaN;

% Salinity is constant
z = -[0 100 500];
s1 = [31 34 35];
sal(:,1) = interp1(z,s1,zprof,'pchip');

% BCstandard
z = -[0 80 200 500];
t1 = [6 0 2 3.5];
temp(:,1) = interp1(z,t1,zprof,'pchip');
% calculate density, for plotting
RHOa = zeros(size(temp)).*NaN;
for i=1:length(temp)
    RHOa(i)=rho(temp(i),sal(i),abs(zprof(i)));
end

% PWcool
t1 = [6 -1 2 3.5];
temp(:,2) = interp1(z,t1,zprof,'pchip');

% PWwarm
t1 = [6 1 2 3.5];
temp(:,3) = interp1(z,t1,zprof,'pchip');

% AWcool
t1 = [6 0 1 2.5];
temp(:,4) = interp1(z,t1,zprof,'pchip');

% AWwarm
t1 = [6 0 3 4.5];
temp(:,5) = interp1(z,t1,zprof,'pchip');

%% Plot initial and open boundary conditions

cmap=parula(4);

% Plot Bathymetry plan view
imagesc(Xdist,Ydist,bathymetry_sill');
colorbar;
xlim([min(Xdist(:)) max(Xdist(:))]);
ylim([min(Ydist(:)) max(Ydist(:))]);
% Model grid
hold on
for i=1:length(Xdist)
    plot([Xdist(i) Xdist(i)],[Ydist(1) Ydist(end)],'k');
end
for i=1:length(Ydist)
    plot([Xdist(1) Xdist(end)],[Ydist(i) Ydist(i)],'k');
end

% Plot bathymetry transect with sill
figure;
plot(Xdist(:),squeeze(bathymetry_sill(:,6)),'color','k','linewidth',2);

% Plot PW conditions
figure;
plot(temp(:,1),zprof(:),'color',[0.7 0.7 0.7],'linewidth',2); hold on; % BCstandard
plot(temp(:,2),zprof(:),'color',cmap(1,:),'linewidth',2); % PWcool
plot(temp(:,3),zprof(:),'color',cmap(3,:),'linewidth',2); % PWwarm
set(gca,'YLim',[-500 0]);
set(gca,'XLim',[-1 6]);

% Plot AW conditions
figure;
plot(temp(:,1),zprof(:),'color',[0.7 0.7 0.7],'linewidth',2); hold on; % BCstandard
plot(temp(:,4),zprof(:),'color',cmap(1,:),'linewidth',2); % AWcool
plot(temp(:,5),zprof(:),'color',cmap(3,:),'linewidth',2); % AWwarm
set(gca,'YLim',[-500 0]);
set(gca,'XLim',[-1 6]);

% Plot standard temperature, salinity and density
figure;
subplot(1,3,1); plot(temp(:,1),zprof,'color',[0.7 0.7 0.7],'linewidth',2);
set(gca,'XLim',[-1 6]);
subplot(1,3,2); plot(sal,zprof,'color',[0.7 0.7 0.7],'linewidth',2);
subplot(1,3,3); plot(RHOa,zprof,'color',[0.7 0.7 0.7],'linewidth',2);



%% Iceberg data

% Concentration
bergConc = zeros(nx,ny,5);
% Scenario 1
conc_linear = linspace(10,1,77); % iceberg concentration declines linearly from 10% adjacent to the glacier to 1 %, over ~40 km.
bergConc(3:length(conc_linear)+2,2:end-1,1) = repmat(conc_linear(:),[1,10]); % iceberg concentration is uniform across fjord
% Scenario 2
conc_linear = linspace(20,1,77); % iceberg concentration declines linearly from 10% adjacent to the glacier to 1 %, over ~40 km.
bergConc(3:length(conc_linear)+2,2:end-1,2) = repmat(conc_linear(:),[1,10]); % iceberg concentration is uniform across fjord
% Scenario 3
conc_linear = linspace(40,1,77); % iceberg concentration declines linearly from 10% adjacent to the glacier to 1 %, over ~40 km.
bergConc(3:length(conc_linear)+2,2:end-1,3) = repmat(conc_linear(:),[1,10]); % iceberg concentration is uniform across fjord
% Scenario 4
conc_linear = linspace(60,5,77); % iceberg concentration declines linearly from 10% adjacent to the glacier to 1 %, over ~40 km.
bergConc(3:length(conc_linear)+2,2:end-1,4) = repmat(conc_linear(:),[1,10]); % iceberg concentration is uniform across fjord
% Scenario 5
conc_linear = linspace(80,5,77); % iceberg concentration declines linearly from 40% adjacent to the glacier to 5 %, over ~40 km.
bergConc(3:length(conc_linear)+2,2:end-1,5) = repmat(conc_linear(:),[1,10]); % iceberg concentration is uniform across fjord

% Maximum iceberg draught
maxDepth = zeros(nx,ny,5);
% Scenario 1
load('iceberg_data/low_discharge_iceberg_data/berg_generation_data','all_berg_depths');
all_depths_1 = all_berg_depths;
% Scenario 2
load('iceberg_data/mid1_discharge_iceberg_data/berg_generation_data','all_berg_depths');
all_depths_2 = all_berg_depths;
% Scenario 3
load('iceberg_data/mid2_discharge_iceberg_data/berg_generation_data','all_berg_depths');
all_depths_3 = all_berg_depths;
% Scenario 4
load('iceberg_data/mid3_discharge_iceberg_data/berg_generation_data','all_berg_depths');
all_depths_4 = all_berg_depths;
% Scenario 5
load('iceberg_data/high_discharge_iceberg_data/berg_generation_data','all_berg_depths');
all_depths_5 = all_berg_depths;
% Map map
for i = 1:nx
    for j = 1:ny
        % Scenario 1
        tmp = all_depths_1{i,j};
        if numel(tmp) > 0
            maxDepth(i,j,1) = nanmax(tmp(:));
        end
        % Scenario 2
        tmp = all_depths_2{i,j};
        if numel(tmp) > 0
            maxDepth(i,j,2) = nanmax(tmp(:));
        end
        % Scenario 3
        tmp = all_depths_3{i,j};
        if numel(tmp) > 0
            maxDepth(i,j,3) = nanmax(tmp(:));
        end
        % Scenario 4
        tmp = all_depths_4{i,j};
        if numel(tmp) > 0
            maxDepth(i,j,4) = nanmax(tmp(:));
        end
        % Scenario 5
        tmp = all_depths_5{i,j};
        if numel(tmp) > 0
            maxDepth(i,j,5) = nanmax(tmp(:));
        end
    end
end

%% Plot icebergs

% Figure
figure('units','normalized','outerposition',[0 0 1 1],'visible','on');  
set(gcf,'color','w');
axes
set(gca, 'color', [1 1 1]);

% Loop
count = 1; % counter for data
for i = 1:10 % loop through plots

    % concentration
    if mod(i,2)~=0
        subplot(5,2,i)
        imagesc(Xdist,Ydist,bergConc(:,:,count)');
        caxis([0 80]);
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
    end
    
    % max draught
    if mod(i,2)==0
        subplot(5,2,i)
        imagesc(Xdist,Ydist,maxDepth(:,:,count)');
        caxis([0 400]);
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        
        % Increment counter
        count = count+1;
    end
    
    % colorbar
    if i == 9
        cb1 = colorbar;
        ylabel(cb1,'Concentration (%)');
        set(gca,'fontsize',12);
    end

    if i == 10
        cb2 = colorbar;
        ylabel(cb2,'Max keel depth (m)');
        set(gca,'fontsize',12);
    end

end

