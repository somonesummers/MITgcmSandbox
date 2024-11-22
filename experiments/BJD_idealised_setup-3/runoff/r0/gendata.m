% This script generates the input files for a MITgcm simulation of an
% idealised fjord, representative of Greenland, utilising the ICEBERG
% package

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


%% Subglacial runoff

% Make masks
runoffVel = zeros(nx,ny,nt);
runoffRad = zeros(nx,ny,nt);
plumeMask = zeros(nx,ny);

% Total runoff (cumecs)
runoff = 0; 

% velocity (m/s) of subglacial runoff
wsg = 1;

% ice front location
icefront=2; % adjacent to wall at western end of domain

% plume location
plume_loc = 6;

%%% Define plume-type mask %%%
% 1 = ice but no plume (melting only)
% 2 = sheet plume (Jenkins)
% 3 = half-conical plume (Morton/Slater)
% 4 = both sheet plume and half-conical plume (NOT YET IMPLEMENTED)
% 5 = detaching conical plume (Goldberg)
% POSITIVE values indicate ice front is orientated north-south
% NEGATIVE values indicate ice front is orientated east-west

% Create virtual ice wall
plumeMask(icefront,2:(end-1)) = 1; % Located 1 cell in from western boundary (need solid barrier behind), and extending across the fjord with (fjord walls either side)
% Specify discharge location
plumeMask(icefront,plume_loc) = 3; % runoff emerges from centre of grounding line

% specify a runoff velocity of 1 m/s
runoffVel(icefront,plume_loc,:) = wsg;

% calculate channel radius
runoffRad(icefront,plume_loc,:) = sqrt(2*runoff/(pi*wsg));

% Write files.
fid=fopen('runoffVel.bin','w','b'); fwrite(fid,runoffVel,acc);fclose(fid);
fid=fopen('runoffRad.bin','w','b'); fwrite(fid,runoffRad,acc);fclose(fid);
fid=fopen('plumeMask.bin','w','b'); fwrite(fid,plumeMask,acc);fclose(fid);


%% Boundary conditions

% pre-allocate
EBCu = zeros(ny,nz,nt);

% Apply barotropic velocity to balance input of runoff
if runoff > 0
    fjordMouthCrossSection = -sum(bathymetry(end,:))*deltaY;
    fjordMouthVelocity = runoff/fjordMouthCrossSection;
    % Out-of-domain velocity is positive at eastern boundary
    EBCu(:) = fjordMouthVelocity;
end

fid=fopen('EBCu.bin','w','b'); fwrite(fid,EBCu,acc);fclose(fid);
