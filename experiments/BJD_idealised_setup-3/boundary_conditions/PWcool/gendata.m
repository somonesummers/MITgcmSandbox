% gendata
% creates boundary conditions based on GLORYS seasonal data (visual
% assesment)

% load data
%load ../../../../../ocean_modification/data/reanalysis/GLORYS/sermilik_mouth/GLORYS_sermilik_seasonal_avg_data.mat


clear
close all
clc

%%  Initial settings

% Add paths to libraries
addpath(genpath('D:/work/PhD/MITgcm/mfiles'));

% Accuracy of binary files
acc = 'real*8';

% Number of time levels for time varying forcing
nt = 1;

% Change in boundary conditions from standard scenario
shift=-1;


%% Gridding

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

% Vertical cell spacing
zprof = -((0.5*deltaZ):deltaZ:((nz*deltaZ)-(0.5*deltaZ)));
delz = zeros(1,nz);
delz(:) = deltaZ;

% Bathymetry
bathymetry = zeros(nx,ny); % pre-allocate
bathymetry(:) = -deltaZ*nz; % uniformly 500 m deep (negative)
bathymetry(1,:) = 0; % barrier behind virtual glacier front (at western end of domain)
bathymetry(:,[1 end]) = 0; % fjord walls


%% Initial conditions

% Below profiles are an idealised version of a profile acquired outside a
% Greenland fjord (based on GLORYS reanalysis data)
z = -[0 80 200 500];
t1 = [6 0+shift 2 3.5];
temp(:,1) = interp1(z,t1,zprof,'pchip');
z = -[0 100 500];
s1 = [31 34 35];
sal(:,1) = interp1(z,s1,zprof,'pchip');

% calculate density, for plotting
RHOa = zeros(size(temp)).*NaN;
for i=1:length(temp)
    RHOa(i)=rho(temp(i),sal(i),abs(zprof(i)));
end

% Plot
subplot(1,3,1);
plot(temp,zprof,'color','r','linewidth',2);
xlabel('Temperature (degC)');
ylabel('Water depth (m)');
title('Initial Conditions');
subplot(1,3,2);
plot(sal,zprof,'color','b','linewidth',2);
xlabel('Salinity (PSU)');
ylabel('Water depth (m)');
subplot(1,3,3);
plot(RHOa,zprof,'color','k','linewidth',2);
xlabel('Density (kg/m3)');
ylabel('Water depth (m)');
print('-dpng','-r150','initial_conditions.png');

% Set initial conditions, uniform throughout domain
saltini = permute(repmat(sal,[1,nx,ny]),[2,3,1]); 
tempini = permute(repmat(temp,[1,nx,ny]),[2,3,1]); 

fid=fopen('saltini.bin','w','b'); fwrite(fid,saltini,acc);fclose(fid);
fid=fopen('tempini.bin','w','b'); fwrite(fid,tempini,acc);fclose(fid);

%% Boundary conditions

% pre-allocate
EBCs = zeros(ny,nz,nt);
EBCt = zeros(ny,nz,nt);

% Make boundary conditions equal to initial conditions
for i = 1:length(temp)
    EBCt(:,i,:) = temp(i);
    EBCs(:,i,:) = sal(i);
end


fid=fopen('EBCs.bin','w','b'); fwrite(fid,EBCs,acc);fclose(fid);
fid=fopen('EBCt.bin','w','b'); fwrite(fid,EBCt,acc);fclose(fid);

