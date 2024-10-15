% This script generates the input files for a MITgcm simulation of an
% idealised fjord, representative of Greenland, utilising the ICEBERG
% package

clear
close all
clc

%%  Initial settings

% Add paths to libraries
addpath(genpath('../../matlabFunctions'));

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
nx=70; % 35 km long
ny=12; % 10 km wide (plus 500 m walls)
nz=20; % 200 m deep

% x scale
delx = zeros(1,nx); 
delx(:) = deltaX;
% fid=fopen('delx.bin','w','b'); fwrite(fid,delx,acc);fclose(fid);

% y scale
dely = zeros(1,ny);
dely(:) = deltaY;
% fid=fopen('dely.bin','w','b'); fwrite(fid,dely,acc);fclose(fid);

% Vertical cell spacing
zprof = -((0.5*deltaZ):deltaZ:((nz*deltaZ)-(0.5*deltaZ)));
delz = zeros(1,nz);
delz(:) = deltaZ;
% 
% Bathymetry
bathymetry = zeros(nx,ny); % pre-allocate
bathymetry(:) = -deltaZ*nz; % uniformly 500 m deep (negative)
% bathymetry(1,:) = 0; % Close western end of domain
bathymetry(:,[1 end]) = 0; % fjord walls
% fid=fopen('bathymetry.bin','w','b'); fwrite(fid,bathymetry,acc);fclose(fid);


%% Initial conditions

% % Below profiles are an idealised version of a profile acquired within a Greenland fjord
% %   - Temperature
% z = -[0 50 200];
% t1 = [0 2 3];
% temp(:,1) = interp1(z,t1,zprof,'linear'); % interpolate to model grid
% %   - Salt
% z = -[0 50 200];
% s1 = [31 34 35];
% sal(:,1) = interp1(z,s1,zprof,'linear'); % interpolate to model grid
% 
% % Set initial conditions, uniform throughout domain
% saltini = permute(repmat(sal,[1,nx,ny]),[2,3,1]); 
% tempini = permute(repmat(temp,[1,nx,ny]),[2,3,1]); 
% fid=fopen('saltini.bin','w','b'); fwrite(fid,saltini,acc);fclose(fid);
% fid=fopen('tempini.bin','w','b'); fwrite(fid,tempini,acc);fclose(fid);
% 
% %% Boundary conditions
% 
% % pre-allocate
% EBCu = zeros(ny,nz,nt);
% EBCs = zeros(ny,nz,nt);
% EBCt = zeros(ny,nz,nt);

% % Make boundary conditions equal to initial conditions
% for i = 1:length(temp)
%     EBCt(:,i,:) = temp(i);
%     EBCs(:,i,:) = sal(i);
% end
% 
% fid=fopen('EBCuBerg.bin','w','b'); fwrite(fid,EBCu,acc);fclose(fid);
% fid=fopen('EBCsBerg.bin','w','b'); fwrite(fid,EBCs,acc);fclose(fid);
% fid=fopen('EBCtBerg.bin','w','b'); fwrite(fid,EBCt,acc);fclose(fid);


%% Icebergs
% Input iceberg files for MITgcm
%   - bergMask: XY grid with +/-1s & 0s for cells containing icebergs. +1
%   where berg long axis is oriented east-west, and -1 where berg long axis
%   is oriented north-south.
%   - bergMaskNums: XY grid with different integers in iceberg cells
%   - numBergsPerCell: XY grid with number of icebergs per cell
%   - driftMask: XY grid of logical 1s & 0s, specifying where to calculate iceberg drift velocity
%   - barrierMask: XY grid of locical 1s & 0s, specifying where to use partial cells to make bergs barriers to water flow
%   - openFrac: XYZ grid specifying fraction of cell that is open (used along with barrierMask)
%   - totalBergArea: XYZ grid specifying total iceberg area in each cell
%   - iceberg_width_n, iceberg_length_n, and iceberg_depth_n : text files
%   containing lists of iceberg width, length and depth in each water column (i.e. cells in plan view. 'n'
%   corresponds to the number in bergMaskNums

% Make masks
bergMask = zeros(nx,ny);
driftMask = zeros(nx,ny);
barrierMask = zeros(nx,ny);
bergVol = zeros(nx,ny);
bergConc = zeros(nx,ny);

% Berg parameters
bergType = 1; % 1 = block; 2 = cone (not implemented)
alpha = 1.8; % slope of inverse power law size frequency distribution
scaling = 2; % 1 = Sulak 2017; 2 = Barker 2004
maxDepth = 75; % (m) - set to zero if 'prescribing' max iceberg width
minDepth= 10; % (m)
maxWidth = 0; % (m) - set to zero if 'prescribing' max iceberg depth
minWidth = 20; % (m)

% Iceberg mask
bergMask(2:31,2:end-1) = 1; % icebergs in inner 5 km, all oriented east-west

% Drift mask
driftMask(2:31,2:end-1) = 1; % calculate effect of iceberg drift on melt rates 

% Barrier mask
barrierMask(2:31,2:end-1) = 1; % make icebergs a physical barrier to water flow

% Iceberg concentration (% of each surface cell that is filled in plan view)
bergConc(2:31,2:end-1) = 75; % iceberg concentration is uniformly 10%

% Generate iceberg size-frequency distribution
 [all_berg_areas, all_berg_lengths, all_berg_widths, all_berg_depths, numBergsPerCell, ...
    total_berg_volume, cell_open_fraction, total_long_face_SA, total_short_face_SA, total_base_SA,...
    total_berg_SA, bergMaskNums ] = ...
    genBerg(bergMask, bergConc, bergType, alpha, scaling, ...
    maxDepth, minDepth, maxWidth, minWidth, ...
    delx, dely, delz, nx, ny, nz, bathymetry);

% write files
fid=fopen('bergMask.bin','w','b'); fwrite(fid,bergMask,acc);fclose(fid);
fid=fopen('bergMaskNums.bin','w','b'); fwrite(fid,bergMaskNums,acc);fclose(fid);
fid=fopen('numBergsPerCell.bin','w','b'); fwrite(fid,numBergsPerCell,acc);fclose(fid);
fid=fopen('openFrac.bin','w','b'); fwrite(fid,cell_open_fraction,acc);fclose(fid);
fid=fopen('driftMask.bin','w','b'); fwrite(fid,driftMask,acc);fclose(fid);
fid=fopen('barrierMask.bin','w','b'); fwrite(fid,barrierMask,acc);fclose(fid);
fid=fopen('totalBergArea.bin','w','b'); fwrite(fid,total_berg_SA,acc);fclose(fid); % m^2

fid=fopen('bergDepths.bin','w','b'); fwrite(fid,total_berg_SA,acc);fclose(fid); % m^2
