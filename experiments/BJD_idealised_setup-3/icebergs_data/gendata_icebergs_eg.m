% gendata_icebergs_eg

clear
close all
clc

%%  Initial settings

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

%% Icebergs
% Final iceberg files for MITgcm
%   - bergMask: XY grid with +/-1s & 0s for cells containing icebergs. +1
%   where berg long axis is oriented east-west, and -1 where berg long axis
%   is oriented north-south.
%   - bergMaskNums: XY grid with different integers in iceberg cells
%   - numBergsPerCell: XY grid with number of icebergs per cell
%   - driftMask: XY grid of logical 1s & 0s, specifying where to calculate iceberg drift velocity
%   - barrierMask: XY grid of locical 1s & 0s, specifying where to use partial cells to make bergs barriers to water flow
%   - openFrac: XYZ grid specifying fraction of cell that is open
%   - totalBergArea: XYZ grid specifying total iceberg area in each cell
%   - iceberg_width_n, iceberg_length_n, and iceberg_depth_n : text files
%   containing lists of iceberg width, length and depth in each cell. 'n'
%   corresponds to the number in bergMaskNums
%
% MITgcm then uses the text files to calculate the surface area of each
% individual iceberg face.

% Make masks
bergMask = zeros(nx,ny);
driftMask = zeros(nx,ny);
barrierMask = zeros(nx,ny);
bergVol = zeros(nx,ny);
bergConc = zeros(nx,ny);

% Berg parameters
bergType = 1; % 1 = block; 2 = cone
alpha = 1.9; % slope of power law size frequency distribution
scaling = 2; % 1 = Sulak 2017; 2 = Barker 2004
maxDepth = 300;
minDepth= 10; 
maxWidth = 0;
minWidth = 20;

% Iceberg mask
bergMask(3:end-20,2:end-1) = 1; % icebergs in every open cell except within 10 km of mouth, all oriented east-west

% Drift mask
driftMask(3:end-20,2:end-1) = 1; % calculate effect of iceberg drift on melt rates 

% Barrier mask
% barrierMask(3:end-20,2:end-1) = 1; % make icebergs a physical barrier to water flow

% Iceberg concentration (% of each surface cell that is filled in plan view)
conc_linear = linspace(80,5,77); % iceberg concentration declines linearly from 40% adjacent to the glacier to 5 %, over ~40 km.
bergConc(3:length(conc_linear)+2,2:end-1) = repmat(conc_linear(:),[1,10]); % iceberg concentration is uniform across fjord

% Generate iceberg size-frequency distribution
 [all_berg_areas, all_berg_lengths, all_berg_widths, all_berg_depths, numBergsPerCell, ...
    total_berg_volume, cell_open_fraction, total_long_face_SA, total_short_face_SA, total_base_SA,...
    total_berg_SA, bergMaskNums ] = ...
    genBerg(bergMask, bergConc, bergType, alpha, scaling, ...
    maxDepth, minDepth, maxWidth, minWidth, ...
    delx, dely, delz, nx, ny, nz, bathymetry);

% show max to check
disp(['The total submerged area is ' num2str(nansum(total_berg_SA(:))/(1000*1000)) ' km squared'])
disp(['The total iceberg volume is ' num2str(nansum(total_berg_volume(:))/(1000*1000*1000)) ' km cubed'])

% write files
fid=fopen('bergMask.bin','w','b'); fwrite(fid,bergMask,acc);fclose(fid);
fid=fopen('bergMaskNums.bin','w','b'); fwrite(fid,bergMaskNums,acc);fclose(fid);
fid=fopen('numBergsPerCell.bin','w','b'); fwrite(fid,numBergsPerCell,acc);fclose(fid);
fid=fopen('openFrac.bin','w','b'); fwrite(fid,cell_open_fraction,acc);fclose(fid);
fid=fopen('driftMask.bin','w','b'); fwrite(fid,driftMask,acc);fclose(fid);
fid=fopen('barrierMask.bin','w','b'); fwrite(fid,barrierMask,acc);fclose(fid);
fid=fopen('totalBergArea.bin','w','b'); fwrite(fid,total_berg_SA,acc);fclose(fid);


% Get variables of interest per cell (plan view)
max_berg_depth_map = zeros(size(bergMask)).*NaN;
min_cell_open_frac_map = zeros(size(bergMask)).*NaN;
num_depths = zeros(size(bergMask)).*NaN;
for i = 1:size(bergMask,1)
    for j = 1:size(bergMask,2)
        if abs(bergMask(i,j)) == 1
            tmp = all_berg_depths{i,j};
            if numel(tmp) > 0
                max_berg_depth_map(i,j) = nanmax(tmp(:));
                num_depths(i,j) = numel(tmp(:));
                tmp = squeeze(cell_open_fraction(i,j,:));
                min_cell_open_frac_map(i,j) = nanmin(tmp(:));
            end
        end
    end
end

% Plot
%   - bergMask
%   - numBergsPerCell
%   - max berg depth per cell
%   - open fraction
subplot(1,4,1); imagesc(bergMask);
subplot(1,4,2); imagesc(numBergsPerCell); cb1 = colorbar; ylabel(cb1,'# bergs');
subplot(1,4,3); imagesc(max_berg_depth_map); cb2 = colorbar; ylabel(cb2, 'max berg depth');
subplot(1,4,4); imagesc(min_cell_open_frac_map); cb3 = colorbar; ylabel(cb3, 'min open frac');


