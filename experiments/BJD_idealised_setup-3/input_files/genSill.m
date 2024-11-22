% genSill

clear
close all
clc

%% Initial settings

% Add paths to libraries
addpath(genpath('D:/work/PhD/MITgcm/mfiles'));

% Accuracy of binary files
acc = 'real*8';

% Number of time levels for time varying forcing
nt = 1;

% sill depth (below fjord surface)
sill_depth = 100; % metres

%% Gridding

% Cell resolution
deltaX = 500;
deltaY = 500;
deltaZ = 10;

% Dimensions of grid
nx=100; % 50 km long
ny=12; % 10 km wide (plus 500 m walls)
nz=50; % 500 m deep

% Sill depth relative to fjord depth
h = sill_depth/(nz*deltaZ);

% x scale
delx = zeros(1,nx); 
delx(:) = deltaX;

% y scale
dely = zeros(1,ny);
dely(:) = deltaY;

% Distances for plotting
Xdist = cumsum(delx./1000);
Ydist = cumsum(dely./1000);

% Vertical cell spacing
zprof = -((0.5*deltaZ):deltaZ:((nz*deltaZ)-(0.5*deltaZ)));
delz = zeros(1,nz);
delz(:) = deltaZ;

% Bathymetry
%   - start with 500 m deep rectangular fjord
bathymetry = zeros(nx,ny); % pre-allocate
bathymetry(:) = -deltaZ*nz; % uniformly 500 m deep (negative)
%   - create gaussian sill
a1 = (deltaZ*nz)-((deltaZ*nz)*h); % amplitude 
b1 = 10; % location of peak along sill (relates to x below)
c1 = 3; % relates to width of peak
x = 1:20; % number of cells wide
sill = a1.*exp(-((x-b1)./c1).^2);
bathymetry(nx-30:nx-30+length(sill)-1,:) = bathymetry(nx-30:nx-30+length(sill)-1,:) + repmat(sill(:),[1,size(bathymetry,2)]);
%   - fjord walls
bathymetry(1,:) = 0; % barrier behind virtual glacier front (at western end of domain)
bathymetry(:,[1 end]) = 0; % fjord walls
fid=fopen('bathymetry.bin','w','b'); fwrite(fid,bathymetry,acc);fclose(fid);

% Plot
figure('units','normalized','outerposition',[0 0 1 1],'visible','on');  
set(gcf,'color','w');
axes
set(gca, 'color', [1 1 1]);
%   - Map
subplot(2,2,1);
imagesc(Xdist,Ydist,bathymetry');
caxis;
cb = colorbar;
ylabel(cb,'Water depth (m)','fontsize',18);
xlabel('Distance along fjord (km)');
ylabel('Distance across fjord (km)');
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
set(gca,'fontsize',18);

% Transect
subplot(2,2,2)
plot(Xdist(:),squeeze(bathymetry(:,6)),'color','k','linewidth',2);
xlabel('Distance along fjord (km)');
ylabel('Depth (m)');
set(gca,'fontsize',18);
export_fig('domain','-r150','-png');
