function [all_berg_areas, all_berg_lengths, all_berg_widths, all_berg_depths, numBergsPerCell, ...
    total_berg_volume, cell_open_fraction, total_long_face_SA, total_short_face_SA, total_base_SA,...
    total_berg_SA, bergMaskNums, bergMask ] =...
    genBerg(bergMask, bergConc, bergType, alpha, scaling, ...
    maxBergDepth, minBergDepth, maxBergWidth, minBergWidth, ...
    delx, dely, delz, nx, ny, nz, bathymetry)
%
% Function to generate iceberg size frequency distributions and iceberg
% submerged surface areas for MITgcm 
%
% It assumes iceberg size-frequency distribution can be described using an inverse power law.
%
% see: https://uk.mathworks.com/matlabcentral/answers/157557-power-law-distribution-how-could-i-generate-random-variables-from-the-truncated-power-law-distribu
% see: https://uk.mathworks.com/matlabcentral/answers/uploaded_files/19033/InverseTransformSampling_InversePowerLaw_RandomDistribution.m
% see: https://uk.mathworks.com/matlabcentral/answers/uploaded_files/19006/rayleigh_random_distribution.m

%% Set up

% Resolution
[X,Y,Z] = meshgrid(delx,dely,delz);
X = permute(X,[2,1,3]);
Y = permute(Y,[2,1,3]);
Z = permute(Z,[2,1,3]);

% Surface area of fjord
cell_SA = X(:,:,1).*Y(:,:,1); % metres squared

% Surface area of fjord covered by icebergs
target_berg_SA = cell_SA.*bergConc./100; % metres squared 

% Calculate average cover
total_target_area = 0;
total_area_of_cells_containing_bergs = 0;
target_area_map = zeros(size(bergMask)).*NaN;
for i = 1:size(bergConc,1)
    for j = 1:size(bergConc,2)
        if abs(bergMask(i,j)) == 1
            total_area_of_cells_containing_bergs = total_area_of_cells_containing_bergs + cell_SA(i,j);
            dummy = (cell_SA(i,j)*(bergConc(i,j)/100));
            total_target_area = total_target_area + dummy;
            target_area_map(i,j) = dummy;
            clear dummy
        end
    end
end
target_percent_cover = total_target_area/total_area_of_cells_containing_bergs*100;

% Depth of water column (at bottom of cell)
depth = cumsum(delz);

% number of cells containing bergs
num_berg_cells = nansum(abs(bergMask(:)));

% plotting
fontSize=12;

% Delete all berg geometry text files
geomList = dir('iceberg_depth*.txt');
for i = 1:numel(geomList)
    delete(geomList(i).name);
end
geomList = dir('iceberg_width*.txt');
for i = 1:numel(geomList)
    delete(geomList(i).name);
end
geomList = dir('iceberg_length*.txt');
for i = 1:numel(geomList)
    delete(geomList(i).name);
end

%% Calculate limits of distribution

if scaling == 1 % then use Sulak17 volume-area scaling
    % volume = 6.0*area^1.30
    if bergType == 1 % block berg
        % assumes volume = L*W*D; and W = L/1.62 (Dowdeswell et al 1992)
        if maxBergWidth==0
            maxBergWidth = 0.0642449*maxBergDepth^(5/3);
            % minBergWidth = 0.0642449*minBergDepth^(5/3);
        elseif maxBergDepth==0
            maxBergDepth = 5.19155*maxBergWidth^(5/3);
            minBergDepth = 5.19155*minBergWidth^(5/3);
        end
    elseif bergType == 2 % cone berg
        if maxBergWidth==0
            maxBergWidth = ((maxBergDepth^(5/3)) / ( 54*2^(2/3)*3^(1/3)*pi^(1/2)))*2;
            minBergWidth = ((minBergDepth^(5/3)) / ( 54*2^(2/3)*3^(1/3)*pi^(1/2)))*2;
        elseif maxBergDepth==0
            maxBergDepth = 9*2^(2/5)*pi^(3/10)*maxBergWidth^(3/5);
            minBergDepth = 9*2^(2/5)*pi^(3/10)*minBergWidth^(3/5);
        end
    end
elseif scaling == 2 % Then use Barker04 width-depth relationship
    % Depth = 2.91*Width^0.71
    if maxBergWidth==0
        maxBergWidth = (100*10^(58/71)*maxBergDepth^(100/71)) / (291*291^(29/71));
        %minBergWidth = (100*10^(58/71)*minBergDepth^(100/71)) / (291*291^(29/71));        
    elseif maxBergDepth==0
        maxBergDepth = 2.91*maxBergWidth^0.71;
        minBergDepth = 2.91*minBergWidth^0.71;
    end
end

%% Generate icebergs
% Create iceberg with dimensions that fit a power law size frequency
% distribution, such that each cell is filled by the desired amount

current_cover = 101; % start with impossibly high cover
number_of_bergs = round(total_target_area / (maxBergWidth*(maxBergWidth/1.62))*30*maxBergDepth/300); % approximate first guess based on ice area, refine later

% Iterate number of bergs and distribution until desired coverage is
% reached
while abs(current_cover - target_percent_cover)>1
    
    clear inversePowerLawDistNumbers_width inversePowerLawDistNumbers_depth
    
    % Iterate number of bergs
    if current_cover > target_percent_cover
        if current_cover-target_percent_cover >= 10
            number_of_bergs = number_of_bergs - round((number_of_bergs*0.05));
        elseif current_cover-target_percent_cover >= 5 && current_cover-target_percent_cover < 10
            number_of_bergs = number_of_bergs - round((number_of_bergs*0.03));
        elseif current_cover-target_percent_cover >= 3 && current_cover-target_percent_cover < 5
            number_of_bergs = number_of_bergs - round((number_of_bergs*0.02));
        else
            number_of_bergs = number_of_bergs - round((number_of_bergs*0.01));
        end
    else
        if target_percent_cover-current_cover >= 10
            number_of_bergs = number_of_bergs + round((number_of_bergs*0.05));
        elseif target_percent_cover-current_cover >= 5 && target_percent_cover-current_cover < 10
            number_of_bergs = number_of_bergs + round((number_of_bergs*0.03));
        elseif target_percent_cover-current_cover >= 3 && target_percent_cover-current_cover < 5
            number_of_bergs = number_of_bergs + round((number_of_bergs*0.02));
        else
            number_of_bergs = number_of_bergs + round((number_of_bergs*0.01));
        end
    end
    
    % Generate the Inverse Power Law cumulative distribution function
    % over the range minBergWidth-maxBergWidth with a slope of alpha.
    x_width = minBergWidth: (maxBergWidth-minBergWidth)/(number_of_bergs*10e2) : maxBergWidth;
    x_depth = minBergDepth: (maxBergDepth-minBergDepth)/(number_of_bergs*10e2) : maxBergDepth;
    inversePowerLawPDF_width = ((alpha-1) / minBergWidth) .* (x_width./minBergWidth) .^ (-alpha);
    inversePowerLawPDF_depth = ((alpha-1) / minBergDepth) .* (x_depth./minBergDepth) .^ (-alpha);
    % Get the CDF numerically
    inversePowerLawCDF_width = cumsum(inversePowerLawPDF_width);
    inversePowerLawCDF_depth = cumsum(inversePowerLawPDF_depth);
    % Normalize
    inversePowerLawCDF_width = inversePowerLawCDF_width ./ inversePowerLawCDF_width(end);
    inversePowerLawCDF_depth = inversePowerLawCDF_depth ./ inversePowerLawCDF_depth(end);
    
    % Generate number_of_bergs uniformly distributed random numbers.
    uniformlyDistributedRandomNumbers = rand(number_of_bergs, 1);
    
    % Invert the CDF of the Inverse Power Law function to get a function that can
    % generate random numbers drawn from a Inverse Power Law distribution,
    % given numbers drawn from a uniform distribution.
    inversePowerLawDistNumbers_width = zeros(1,length(uniformlyDistributedRandomNumbers));
    inversePowerLawDistNumbers_depth = zeros(1,length(uniformlyDistributedRandomNumbers));
    [~,nearestIndex_width] = ismembertol(uniformlyDistributedRandomNumbers(:),inversePowerLawCDF_width(:),nanmax(diff(inversePowerLawCDF_width(:))));
    [~,nearestIndex_depth] = ismembertol(uniformlyDistributedRandomNumbers(:),inversePowerLawCDF_depth(:),nanmax(diff(inversePowerLawCDF_depth(:))));
    inversePowerLawDistNumbers_width = x_width(nearestIndex_width);
    inversePowerLawDistNumbers_depth = x_depth(nearestIndex_depth); 
%     for k = 1:length(uniformlyDistributedRandomNumbers)
%         nearestIndex_width = find(inversePowerLawCDF_width(:) >= uniformlyDistributedRandomNumbers(:), 1, 'first');
%         nearestIndex_depth = find(inversePowerLawCDF_depth >= uniformlyDistributedRandomNumbers(k), 1, 'first');
%         inversePowerLawDistNumbers_width(k) = x_width(nearestIndex_width);
%         inversePowerLawDistNumbers_depth(k) = x_depth(nearestIndex_depth);    
%     end

    % Check percent coverage
    bergArea = inversePowerLawDistNumbers_width.*(inversePowerLawDistNumbers_width./1.62);
    bergArea = sum(bergArea(:));
    current_cover = (bergArea/total_area_of_cells_containing_bergs)*100;
    
    % display current iteration output
    disp(['Target coverage: ' num2str(target_percent_cover) ]);
    disp(['Current coverage: ' num2str(current_cover) ]);
    disp(['Number of icebergs: ' num2str(number_of_bergs) ]);
    
end

%% Plot
% Plot the PDF.
figure('units','normalized','outerposition',[0 0 1 1],'visible','on');  
subplot(3, 2, 1);
plot(x_depth, inversePowerLawPDF_depth, 'b-', 'LineWidth', 3);
caption = sprintf('Depth. Inverse Power Law PDF with alpha = %.2f', alpha);
grid on;
title(caption, 'FontSize', fontSize);
%set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
set(gcf,'name','Iceberg depth-frequency distribution','numbertitle','off') 
% Plot the CDF.
% 
subplot(3, 2, 2);
plot(x_depth, inversePowerLawCDF_depth, 'b-', 'LineWidth', 3);
caption = sprintf('Depth. Inverse Power Law CDF with alpha = %.2f', alpha);
grid on;
title(caption, 'FontSize', fontSize);

subplot(3,2,3);
bar(uniformlyDistributedRandomNumbers, 'BarWidth', 1);
xlim([0 number_of_bergs]);
caption = sprintf('%d Uniformly Distributed Numbers', number_of_bergs);
title(caption, 'FontSize', fontSize);
xlabel('Element Number');
ylabel('Value');
% 
% Plot the Inverse Power Law distributed numbers.
subplot(3,2,4);
bar(inversePowerLawDistNumbers_depth, 'BarWidth', 1);
xlim([0 number_of_bergs]);
caption = sprintf('Depth. %d Inverse Power Law Distributed Numbers', number_of_bergs);
title(caption, 'FontSize', fontSize);
xlabel('Element Number');
ylabel('Value');
% 
% Get histogram of uniformly distributed numbers.
[countsU, binsU] = hist(uniformlyDistributedRandomNumbers, 50);
% Plot the uniformly distributed numbers.
subplot(3,2,5);
bar(binsU, countsU, 'BarWidth', 1);
grid on;
caption = sprintf('Histogram of %d Uniformly Distributed Numbers', number_of_bergs);
title(caption, 'FontSize', fontSize);
xlabel('Value');
ylabel('Count');
% 
% Get histogram of Inverse Power Law distributed numbers.
% Observe that it's distribution is not flat like it is
% for the uniformly distributed numbers.
% It will take on the Inverse Power Law distribution shape.
[countsR, binsR] = hist(inversePowerLawDistNumbers_depth, 50);
subplot(3,2,6);
bar(binsR, countsR, 'BarWidth', 1);
grid on;
%caption = sprintf('Width. Histogram of %d Inverse Power Law Distributed Numbers', number_of_bergs);
caption = 'Histogram of bergDepth';
title(caption, 'FontSize', fontSize);
xlabel('Value');
ylabel('Count');
print('-dpng','-r150','berg_depth_frequency_distribution.png');
close

%% Plot width
% Plot the PDF.
figure('units','normalized','outerposition',[0 0 1 1],'visible','on');  
subplot(3, 2, 1);
plot(x_depth, inversePowerLawPDF_width, 'b-', 'LineWidth', 3);
caption = sprintf('Width. Inverse Power Law PDF with alpha = %.2f', alpha);
grid on;
title(caption, 'FontSize', fontSize);
%set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
set(gcf,'name','Iceberg length-frequency distribution','numbertitle','off') 
% Plot the CDF.
% 
subplot(3, 2, 2);
plot(x_depth, inversePowerLawCDF_width, 'b-', 'LineWidth', 3);
caption = sprintf('Width. Inverse Power Law CDF with alpha = %.2f', alpha);
grid on;
title(caption, 'FontSize', fontSize);

subplot(3,2,3);
bar(uniformlyDistributedRandomNumbers, 'BarWidth', 1);
xlim([0 number_of_bergs]);
caption = sprintf('%d Uniformly Distributed Numbers', number_of_bergs);
title(caption, 'FontSize', fontSize);
xlabel('Element Number');
ylabel('Value');
% 
% Plot the Inverse Power Law distributed numbers.
subplot(3,2,4);
bar(inversePowerLawDistNumbers_width, 'BarWidth', 1);
xlim([0 number_of_bergs]);
caption = sprintf('Width. %d Inverse Power Law Distributed Numbers', number_of_bergs);
title(caption, 'FontSize', fontSize);
xlabel('Element Number');
ylabel('Value');
% 
% Get histogram of uniformly distributed numbers.
[countsU, binsU] = hist(uniformlyDistributedRandomNumbers, 50);
% Plot the uniformly distributed numbers.
subplot(3,2,5);
bar(binsU, countsU, 'BarWidth', 1);
grid on;
caption = sprintf('Histogram of %d Uniformly Distributed Numbers', number_of_bergs);
title(caption, 'FontSize', fontSize);
xlabel('Value');
ylabel('Count');
% 
% Get histogram of Inverse Power Law distributed numbers.
% Observe that it's distribution is not flat like it is
% for the uniformly distributed numbers.
% It will take on the Inverse Power Law distribution shape.
[countsR, binsR] = hist(inversePowerLawDistNumbers_width, 50);
subplot(3,2,6);
bar(binsR, countsR, 'BarWidth', 1);
grid on;
%caption = sprintf('Width. Histogram of %d Inverse Power Law Distributed Numbers', number_of_bergs);
caption = 'Histogram of bergWidth';
title(caption, 'FontSize', fontSize);
xlabel('Value');
ylabel('Count');
print('-dpng','-r150','berg_width_frequency_distribution.png');
close

%% Calculate depth-invariant aspects of iceberg geometry
% Calculate geometry of entire berg
if bergType == 1
    bergLongWidth = inversePowerLawDistNumbers_width;
    bergShortWidth(:) = bergLongWidth(:)./1.62; 
    bergDepth = inversePowerLawDistNumbers_depth;
    bergTopArea = bergLongWidth(:).*bergShortWidth(:);
elseif bergType == 2
    bergTopRadius = inversePowerLawDistNumbers_width(:)./2;
    bergTopArea = pi.*bergTopRadius(:).^2;
    bergDepth = inversePowerLawDistNumbers_depth;
   % bergSlLength = sqrt((bergTopRadius(:).^2)+(inversePowerLawDistNumbers_depth(:).^2));
    theta = 90 - (atand(inversePowerLawDistNumbers_depth(:)./bergTopRadius(:)));
end

% Some variables are constant through the water column
if bergType == 2
    tslant = zeros(nz,length(bergLongWidth)); % length of sloping side of cone in a single cell
    bergLRadius = zeros(nz,length(bergLongWidth)); % Radius of cone at base of cell
    tslant_ice = 0; % iceTopoConst./cosd(theta(:)); % length of sloping side of cone in sea ice
    bergLRadius_ice = bergTopRadius(:) - (tslant_ice.*sind(theta(:)));  % length of sloping side of cone in sea ice
    % now loop through water column, in case vertical resolution varies
    for i = 1:nz
        tslant(i,:) = delz(i)./cosd(theta(:));
        if i == 1
            bergLRadius(i,:) = bergTopRadius(:)' - (tslant(i,:).*sind(theta(:)'));
            bergLRadius(bergLRadius<0) = 0;
        else
            bergLRadius(i,:) = bergLRadius(i-1,:) - (tslant(i,:).*sind(theta(:)'));
            bergLRadius(bergLRadius<0) = 0;
        end
    end
end
clear i

%% Deal with large bergs

% 1. First, loop through cells and identify large bergs
all_large_bergs_index = cell(size(bergMask));
count = 1;
for i = 1:size(bergMask,1)
    for j = 1:size(bergMask,2)
        
        % only bother if there are bergs in this cell
        if abs(bergMask(i,j)) == 1
            
            % Surface area of current cell
            local_cell_SA = cell_SA(i,j);
            
            % Identify bergs with at least one axis longer than this cell
            %l = a > local_cell_SA;
            if bergMask(i,j) == 1
                l = bergLongWidth > X(i,j,1); % assumes X is oriented east-west
            elseif bergMask(i,j) == -1
                l = bergLongWidth > Y(i,j,1); % assumes Y is oriented north-south
            end
            all_large_bergs_index{i,j} = find(l==1);
            
            % Also keep track of linear index of cells that could contain a
            % large berg
            if nansum(l(:)) > 0
                potential_cell_index(count) = sub2ind(size(bergMask),i,j);
                count = count+1;
            end
            
        end
        
    end
end

% Make list of bergs that need to be allocated
large_berg_list = unique(horzcat(all_large_bergs_index{:})); % indexes that point to the large bergs


% 2. Next, send identified large bergs to specific cells
% Get linear indices of cells that could contain large icebergs
if exist('potential_cell_index','var')
    r = potential_cell_index(randperm(length(potential_cell_index)));

    % Allocate large bergs to cells
    allocated_bergs = zeros(size(large_berg_list)).*NaN;
    location_of_large_bergs = zeros(size(large_berg_list)).*NaN;
    all_allocated_bergs = cell(size(bergMask));
    count = 1;
    allocation_count = 1;
    while nansum(~isnan(allocated_bergs(:))) < numel(large_berg_list(:))

        % Get row/column of random cell
        [row,col] =  ind2sub(size(bergMask),r(count));

        % Get list of large bergs in this cell
        local_large_bergs = all_large_bergs_index{row,col};

        % Allocate one of the unused large bergs to this cell, if one is
        % available
        avail = ismember(local_large_bergs,allocated_bergs); % logical 0 where the berg has not already been allocated
        % Only continue if there is at least one berg that could go in this cell
        if ismember(0,avail)
            berg_to_use = local_large_bergs(find(avail == 0, 1, 'first')); % get index of berg to allocate
            tmp = all_allocated_bergs{row,col}; % get list of bergs already in this cell (there shouldn't be any, but just in case)
            all_allocated_bergs{row,col} = [tmp(:);berg_to_use(:)]; % append large bergs to this cell
            allocated_bergs(allocation_count) = berg_to_use; % keep track of allocated bergs
            location_of_large_bergs(allocation_count) = sub2ind(size(bergMask),row,col); % keep track of linear index of berg locations
            allocation_count = allocation_count+1; % increment counter
        end

        % Increment counter
        count = count+1;

    end


    % 3. Split up large bergs
    for i = 1:length(location_of_large_bergs)

        % Get cell location
        ind = location_of_large_bergs(i);
        [row,col] = ind2sub(size(bergMask),ind);

        % Get berg index
        berg_ind = allocated_bergs(i);

        % Get berg dimensions
        l = bergLongWidth(berg_ind);
        w = bergShortWidth(berg_ind);
        d = bergDepth(berg_ind);

        % Identify which dimension is cell long axis
        if bergMask(row,col) == 1
            long_dim = X;
            short_dim = Y;
        else
            long_dim = Y;
            short_dim = X;
        end

        % How many cells wide and long is the iceberg from the berg centre
        % (assumes constant resolution!)
        %   -1 = berg does not extent out of the central cell in that direction
        %   0 = berg extends out of the central cell in that direction
        %   n = berg extends out of the central cell and n other cells in that direction
        row_buffer = round((l/2 - long_dim(row,col,1))/long_dim(row,col,1));
        col_buffer = round((w/2 - short_dim(row,col,1))/short_dim(row,col,1));

        % Set range of cells to consider (to simplify code)
        row_range = row-max([0,row_buffer]):row+max([0,row_buffer]);
        col_range = col-max([0,col_buffer]):col+max([0,col_buffer]);

        % Get dimensions of interest
        %   - central cell (or cells not containing berg edges)
        width_C = nanmin( repmat(w,[length(row_range),length(col_range)]) , short_dim(row_range , col_range,1)); % width in centre cell (could be iceberg or cell width)
        length_C = long_dim(row_range , col_range,1);  % length of cell long axis (cell length will always be less than berg length, so no need to compare
        %   - northern berg edge cell(s)
        width_N = nanmin( repmat(w,[length(col_range),1]) , short_dim(min(row_range) , col_range,1)'); % width in northern cell is the same as central cell (could be iceberg or cell width)
        length_N = repmat((l - sum(long_dim(row_range,col))) / 2 , [length(col_range),1]); % length in northern cell: half of the berg length minus the length of contained cells
        %   - southern berg edge cell(s)
        width_S = nanmin( repmat(w,[length(col_range),1]) , short_dim(min(row_range) , col_range,1)'); % width in southern cell is the same as central cell (could be iceberg or cell width)
        length_S = repmat((l - sum(long_dim(row_range,col))) / 2 , [length(col_range),1]); % length in southern cell: half of the berg length minus the length of contained cells
        %   - eastern berg edge cell(s)
        width_E = repmat(nanmax([0, (w - sum(short_dim(row,col_range))) / 2]), [length(row_range),1]); % width in eastern cell: 0 if the berg is not wider than the cell, or half the difference between the contained cell widths and the berg width
        length_E = long_dim(row_range,max(col_range)); % length of cell long axis
        %   - western cell
        width_W = repmat(nanmax([0, (w - sum(short_dim(row,col_range))) / 2]), [length(row_range),1]); % width in western cell: 0 if the berg is not wider than the cell, or half the difference between the contained cell widths and the berg width
        length_W = long_dim(row_range,min(col_range)); % length of cell long axis
        %   - northwestern cell
        width_NW = nanmax([0, (w - sum(short_dim(row,col_range))) / 2]); % width in northwestern cell: 0 if the berg is not wider than the cell, or half the difference between the contained cell widths and the berg width
        length_NW = (l - sum(long_dim(row_range,col))) / 2; % length in northwestern cell: half of the berg length minus the length of contained cells
        %   - northeastern cell
        width_NE = nanmax([0, (w - sum(short_dim(row,col_range))) / 2]); % width in northeastern cell: 0 if the berg is not wider than the cell, or half the difference between the contained cell widths and the berg width
        length_NE = (l - sum(long_dim(row_range,col))) / 2; % length in northeastern cell: half of the berg length minus the length of contained cells
        %   - southeastern cell
        width_SE = nanmax([0, (w - sum(short_dim(row,col_range))) / 2]); % width in southeastern cell: 0 if the berg is not wider than the cell, or half the difference between the contained cell widths and the berg width
        length_SE = (l - sum(long_dim(row_range,col))) / 2; % length in southeastern cell: half of the berg length minus the length of contained cells
        %   - southwestern cell
        width_SW = nanmax([0, (w - sum(short_dim(row,col_range))) / 2]); % width in southwestern cell: 0 if the berg is not wider than the cell, or half the difference between the contained cell widths and the berg width
        length_SW = (l - sum(long_dim(row_range,col))) / 2; % length in southwestern cell: half of the berg length minus the length of contained cells

        % Now split these up into separate icebergs? 
        %   - central cells
        bergLongWidth(end+1:end+length(length_C)) = length_C;
        bergShortWidth(end+1:end+length(width_C)) = width_C;
        bergTopArea(end+1:end+length(width_C)) = length_C.*width_C;
        bergDepth(end+1:end+length(width_C)) = d;
        %   - north/south cells
        bergLongWidth(end+1:end+length(length_N)) = length_N;
        bergShortWidth(end+1:end+length(width_N)) = width_N;
        bergTopArea(end+1:end+length(width_N)) = length_N.*width_N;
        bergDepth(end+1:end+length(width_N)) = d;
        bergLongWidth(end+1:end+length(length_S)) = length_S;
        bergShortWidth(end+1:end+length(width_S)) = width_S;
        bergTopArea(end+1:end+length(width_S)) = length_S.*width_S;
        bergDepth(end+1:end+length(width_S)) = d;
        %   - east/west cells
        count = 1;
        for r = row_range
            if width_E(count) > 0 && length_E(count) > 0 % only bother if cell actually extends into this cell
                bergLongWidth(end+1) = length_E(count);
                bergShortWidth(end+1) = width_E(count);
                bergTopArea(end+1) = length_E(count)*width_E(count);
                bergDepth(end+1) = d;
                bergLongWidth(end+1) = length_W(count);
                bergShortWidth(end+1) = width_W(count);
                bergTopArea(end+1) = length_W(count)*width_W(count);
                bergDepth(end+1) = d;
                count = count+1;
            end
        end
        %   - northeast cell
        if width_NE > 0 && length_NE > 0 % only bother if both dimensions are greater than zero
            bergLongWidth(end+1) = length_NE;
            bergShortWidth(end+1) = width_NE;
            bergTopArea(end+1) = length_NE*width_NE;
            bergDepth(end+1) = d;
        end
        %   - northwest cell
        if width_NW > 0 && length_NW > 0 % only bother if both dimensions are greater than zero
            bergLongWidth(end+1) = length_NW;
            bergShortWidth(end+1) = width_NW;
            bergTopArea(end+1) = length_NW*width_NW;
            bergDepth(end+1) = d;
        end
        %   - southeast cell
        if width_SE > 0 && length_SE > 0 % only bother if both dimensions are greater than zero
            bergLongWidth(end+1) = length_SE;
            bergShortWidth(end+1) = width_SE;
            bergTopArea(end+1) = length_SE*width_SE;
            bergDepth(end+1) = d;
        end
        %   - southwest cell
        if width_SW > 0 && length_SW > 0 % only bother if both dimensions are greater than zero
            bergLongWidth(end+1) = length_SW;
            bergShortWidth(end+1) = width_SW;
            bergTopArea(end+1) = length_SW*width_SW;
            bergDepth(end+1) = d;
        end

    end

    % Remove large bergs from original list
    bergLongWidth(large_berg_list) = [];
    bergShortWidth(large_berg_list) = [];
    bergTopArea(large_berg_list) = [];
    bergDepth(large_berg_list) = [];

    disp(['Number of icebergs after breaking up large bergs: ' num2str(length(bergTopArea)) ]);
    
end

%% Send all bergs to cells, such that cells are appropriately filled

% Setup
a = bergTopArea(:); 
w = bergShortWidth(:);
l = bergLongWidth(:);
d = bergDepth(:);

% pre-allocate
all_berg_areas = cell(size(bergMask));
total_berg_area = zeros(size(bergMask)).*NaN;
all_berg_lengths = cell(size(bergMask));
all_berg_widths = cell(size(bergMask));
all_berg_depths = cell(size(bergMask));
all_allocated_bergs = cell(size(bergMask));
bergMaskNums = zeros(size(bergMask));

% Initiate berg cell counter
berg_cell_count = 1; % This allocates a number to each column containing bergs

% Loop through berg cells randomly (seems to result in closer reproduction
% of target coverage than looping through in one direction)
bm = randperm(length(bergMask(:)));
disp('Assigning bergs to cells...')
for r = 1:length(bm)
    
    [i,j] = ind2sub(size(bergMask),bm(r));
        
    % only bother if there are bergs in this cell and there are still
    % bergs available for allocation
    if abs(bergMask(i,j)) == 1
            
        % Assign cell a number
        bergMaskNums(i,j) = berg_cell_count;
        berg_cell_count = berg_cell_count + 1;

       % Get string for berg dimension text file
       if bergMaskNums(i,j) < 10
           bergNumStr = ['0000' num2str(bergMaskNums(i,j)) ];
       elseif bergMaskNums(i,j) >= 10 && bergMaskNums(i,j) < 100
           bergNumStr = ['000' num2str(bergMaskNums(i,j)) ];
       elseif bergMaskNums(i,j) >= 100 && bergMaskNums(i,j) < 1000
           bergNumStr = ['00' num2str(bergMaskNums(i,j)) ]; 
       elseif bergMaskNums(i,j) >= 1000 && bergMaskNums(i,j) < 10000
           bergNumStr = ['0' num2str(bergMaskNums(i,j)) ];
       else
           bergNumStr = num2str(bergMaskNums(i,j));
       end        
        
        % Make copy of berg dimensions
        a_copy = a;   
        w_copy = w;
        l_copy = l;
        d_copy = d;

        % Keep only non-nan elements, and keep track of locations
        non_nan_idx = find(~isnan(a_copy(:)));
        a_use = a_copy(non_nan_idx(:));
        w_use = w_copy(non_nan_idx(:));
        l_use = l_copy(non_nan_idx(:));
        d_use = d_copy(non_nan_idx(:));

        % Surface area of current cell
        local_cell_SA = cell_SA(i,j);

        % Target % coverage of icebergs in this cell
        local_target_berg_conc = bergConc(i,j);

        % Target surface area of icebergs in this cell
        local_target_berg_SA = local_cell_SA*(local_target_berg_conc/100);

        % How much area do we still need
        remainder = local_target_berg_SA;

        % while the remainder is greater than 1 % of the berg surface area,
        % add or remove bergs as necessary 
        exit_while = 0;
        local_berg_areas = [];
        local_berg_widths = [];
        local_berg_lengths = [];
        local_berg_depths = [];
        while abs(remainder) > local_cell_SA*0.01 && exit_while == 0

            % Find single berg that is close to filling the space
            too_deep = find(abs(d_use) > abs(bathymetry(i,j)));
            a_tmp = a_use;
            a_tmp(too_deep) = NaN;
            difference = abs(a_tmp - remainder);
            clear a_tmp
            idx2 = find(difference == nanmin(difference),1,'first'); 

            % Only append berg if it makes the remainder smaller
            tmp_berg_area = [local_berg_areas(:);a_use(idx2)]; %  Create temporary array to check if new berg helps
            total_area_of_tmp_bergs = nansum(tmp_berg_area(:)); % Get area of bergs with new addition
            new_remainder = local_target_berg_SA - total_area_of_tmp_bergs; % Get new remainder
            if abs(new_remainder) < abs(remainder)
                % append new berg
                local_berg_areas = [local_berg_areas(:);a_use(idx2)];
                local_berg_widths = [local_berg_widths(:);w_use(idx2)];
                local_berg_lengths = [local_berg_lengths(:);l_use(idx2)];
                local_berg_depths = [local_berg_depths(:);d_use(idx2)];
                % set new berg to NaN, so it is not repeated
                % - master copy
                used_berg = find(a_copy(:) == a_use(idx2));
                a_copy(used_berg) = NaN;
                w_copy(used_berg) = NaN;
                l_copy(used_berg) = NaN;
                d_copy(used_berg) = NaN;
                %   - working copy
                a_use(idx2) = NaN; 
                w_use(idx2) = NaN;
                l_use(idx2) = NaN;
                d_use(idx2) = NaN;    
                % assign remainder
                remainder = new_remainder;
                % clear variables
                clear difference idx2 tmp_berg_area total_area_of_tmp_bergs new_remainder
            else
                exit_while = 1;
            end
            clear tmp_berg_area total_area_of_tmp_bergs new_remainder idx2 difference

        end

        % Save bergs to cell array with location
        tmp = all_allocated_bergs{i,j};
        [~,tmp2] = intersect(bergTopArea(:), local_berg_areas(:),'stable');
        all_allocated_bergs{i,j} = [tmp(:) ; tmp2(:)];
        tmp = all_berg_areas{i,j};
        all_berg_areas{i,j} = [tmp(:);local_berg_areas(:)]; % list of area of individual bergs in each cell
        tmp = all_berg_widths{i,j};
        all_berg_widths{i,j} = [tmp(:);local_berg_widths(:)];
        tmp = all_berg_lengths{i,j};
        all_berg_lengths{i,j} = [tmp(:);local_berg_lengths(:)];
        tmp = all_berg_depths{i,j};
        all_berg_depths{i,j} = [tmp(:);local_berg_depths(:)];
        clear tmp

        % Convert area to double
        dummy = all_berg_areas{i,j};
        total_berg_area(i,j) = nansum(dummy(:));

        % Overwrite 'a' so that bergs that have been used here aren't
        % repeated elsewhere
        a = a_copy;
        w = w_copy;
        l = l_copy;
        d = d_copy;
        
        % Now that this cell is done, write the lists to a file
        %   - lists of iceberg depth in each cell
        %   - lists of iceberg width in each cell
        %   - lists of iceberg length in each cell
        % Files named e.g. iceberg_depth_1, iceberg_depth_2, ... iceberg_depth_n
        % where the number relates to the number assigned to that cell (provided by
        % bergMaskNums)
        % delete existing files with same name
        if exist(['iceberg_width_' bergNumStr '.txt'],'file') == 2
            delete(['iceberg_width_' bergNumStr '.txt']);
            delete(['iceberg_length_' bergNumStr '.txt']);
            delete(['iceberg_depth_' bergNumStr '.txt']);
        end
        % get list
        width_list = all_berg_widths{i,j};
        depth_list = all_berg_depths{i,j};
        length_list = all_berg_lengths{i,j};
        if numel(width_list) > 0
            % open files
            fid_width = fopen(['iceberg_width_' bergNumStr '.txt'],'w');
            fid_length = fopen(['iceberg_length_' bergNumStr '.txt'],'w');
            fid_depth = fopen(['iceberg_depth_' bergNumStr '.txt'],'w');
            % write to files
            fprintf(fid_width,'%f\n',width_list);
            fprintf(fid_length,'%f\n',length_list);
            fprintf(fid_depth,'%f\n',depth_list);
            % close files
            fclose(fid_width); 
            fclose(fid_length);
            fclose(fid_depth);
            % clear
            clear fid_width fid_length fid_depth
        end
            
    end % check berg mask
        
end % loop through cells

% See how we did
derived_cover = total_berg_area./cell_SA.*100;
difference = derived_cover - bergConc;

% Get cells where there is space for more bergs
idx = find(difference<0);

% Loop through cells again and use small icebergs? 
if numel(idx) > 0 && sum(~isnan(a_copy)) > 0
    disp('Assigning small bergs to cells...')
    used_all_bergs = 0;
    for r = 1:length(idx)

        % Get row and column index of this cell
        [i,j] = ind2sub(size(bergMask),idx(r));

        % Get derived and target berg areas
        derived_area = nansum(all_berg_areas{i,j});
        target_area = cell_SA(i,j)*(bergConc(i,j)/100);

        % only bother if there are bergs in this cell and there are still
        % bergs available for allocation
        if abs(bergMask(i,j)) == 1 && derived_area < target_area && used_all_bergs == 0
            
            % Get string for berg dimension text file
           if bergMaskNums(i,j) < 10
               bergNumStr = ['0000' num2str(bergMaskNums(i,j)) ];
           elseif bergMaskNums(i,j) >= 10 && bergMaskNums(i,j) < 100
               bergNumStr = ['000' num2str(bergMaskNums(i,j)) ];
           elseif bergMaskNums(i,j) >= 100 && bergMaskNums(i,j) < 1000
               bergNumStr = ['00' num2str(bergMaskNums(i,j)) ]; 
           elseif bergMaskNums(i,j) >= 1000 && bergMaskNums(i,j) < 10000
               bergNumStr = ['0' num2str(bergMaskNums(i,j)) ];
           else
               bergNumStr = num2str(bergMaskNums(i,j));
           end    

            % How much area is left to fill?
            remainder = target_area - derived_area;

            % set-up
            exit_while = 0;
            local_berg_areas = [];
            local_berg_widths = [];
            local_berg_lengths = [];
            local_berg_depths = [];

            % Make copy of berg dimensions
            a_copy = a;   
            w_copy = w;
            l_copy = l;
            d_copy = d;

            % Keep only non-nan elements, and keep track of locations
            non_nan_idx = find(~isnan(a_copy(:)));
            if numel(non_nan_idx) == 0
                used_all_bergs = 1;
                break
            end
            a_use = a_copy(non_nan_idx(:));
            w_use = w_copy(non_nan_idx(:));
            l_use = l_copy(non_nan_idx(:));
            d_use = d_copy(non_nan_idx(:));

            % Add small bergs until the area is full or until all bergs are
            % used?
             while exit_while == 0 && used_all_bergs == 0

                % Find the smallest berg that will fit in the space?
                too_deep = find(abs(d_use) > abs(bathymetry(i,j)));
                a_tmp = a_use;
                a_tmp(too_deep) = NaN;
                difference = (a_tmp - remainder).*-1; % positive where the berg is smaller than the remaining space
                if nanmax(difference(:)) > 0
                    idx2 = find(difference == nanmax(difference),1,'first'); 
                else
                    exit_while = 1;
                    difference = abs(a_tmp - remainder);
                    idx2 = find(difference == nanmin(difference),1,'first'); 
                end
                
                if numel(idx2) > 0

                % Only append berg if it makes the remainder smaller
                tmp_berg_area = [local_berg_areas(:);a_use(idx2)]; %  Create temporary array to check if new berg helps
                total_area_of_tmp_bergs = nansum(tmp_berg_area(:)); % Get area of bergs with new addition
                new_remainder = target_area - total_area_of_tmp_bergs; % Get new remainder
                if new_remainder < remainder || exit_while == 1
                    % append new berg
                    local_berg_areas = [local_berg_areas(:);a_use(idx2)];
                    local_berg_widths = [local_berg_widths(:);w_use(idx2)];
                    local_berg_lengths = [local_berg_lengths(:);l_use(idx2)];
                    local_berg_depths = [local_berg_depths(:);d_use(idx2)];
                    % set new berg to NaN, so it is not repeated
                    % - master copy
                    used_berg = find(a_copy(:) == a_use(idx2));
                    a_copy(used_berg) = NaN;
                    w_copy(used_berg) = NaN;
                    l_copy(used_berg) = NaN;
                    d_copy(used_berg) = NaN;
                    %   - working copy
                    a_use(idx2) = NaN; 
                    w_use(idx2) = NaN;
                    l_use(idx2) = NaN;
                    d_use(idx2) = NaN;    
                    % assign remainder
                    remainder = new_remainder;
                    % clear variables
                    clear difference idx2 tmp_berg_area total_area_of_tmp_bergs new_remainder
                else
                    exit_while = 1;
                end
                clear tmp_berg_area total_area_of_tmp_bergs new_remainder idx2 difference
                
                end

            end

            % Save bergs to cell array with location
            tmp = all_allocated_bergs{i,j};
            [~,tmp2] = intersect(bergTopArea(:), local_berg_areas(:),'stable');
            all_allocated_bergs{i,j} = [tmp(:) ; tmp2(:)];
            tmp = all_berg_areas{i,j};
            all_berg_areas{i,j} = [tmp(:);local_berg_areas(:)]; % list of area of individual bergs in each cell
            tmp = all_berg_widths{i,j};
            all_berg_widths{i,j} = [tmp(:);local_berg_widths(:)];
            tmp = all_berg_lengths{i,j};
            all_berg_lengths{i,j} = [tmp(:);local_berg_lengths(:)];
            tmp = all_berg_depths{i,j};
            all_berg_depths{i,j} = [tmp(:);local_berg_depths(:)];
            clear tmp

            % Convert area to double
            dummy = all_berg_areas{i,j};
            total_berg_area(i,j) = nansum(dummy(:));

            % Overwrite 'a' so that bergs that have been used here aren't
            % repeated elsewhere
            a = a_copy;
            w = w_copy;
            l = l_copy;
            d = d_copy;

            % Now that this cell is done, write the lists to a file
            %   - lists of iceberg depth in each cell
            %   - lists of iceberg width in each cell
            %   - lists of iceberg length in each cell
            % Files named e.g. iceberg_depth_1, iceberg_depth_2, ... iceberg_depth_n
            % where the number relates to the number assigned to that cell (provided by
            % bergMaskNums)
            % delete existing files with same name
            if exist(['iceberg_width_' bergNumStr '.txt'],'file') == 2
                delete(['iceberg_width_' bergNumStr '.txt']);
                delete(['iceberg_length_' bergNumStr '.txt']);
                delete(['iceberg_depth_' bergNumStr '.txt']);
            end
            % get list
            width_list = all_berg_widths{i,j};
            depth_list = all_berg_depths{i,j};
            length_list = all_berg_lengths{i,j};
            % open files
            if numel(width_list) > 0
                fid_width = fopen(['iceberg_width_' bergNumStr '.txt'],'w');
                fid_length = fopen(['iceberg_length_' bergNumStr '.txt'],'w');
                fid_depth = fopen(['iceberg_depth_' bergNumStr '.txt'],'w');
                % write to files
                fprintf(fid_width,'%f\n',width_list);
                fprintf(fid_length,'%f\n',length_list);
                fprintf(fid_depth,'%f\n',depth_list);
                % close files
                fclose(fid_width); 
                fclose(fid_length);
                fclose(fid_depth);
                % clear
                clear fid_width fid_length fid_depth
            end

        end % check berg mask

    end % loop through cells
end % check if there is space for more bergs

% Plot output to check
derived_cover = total_berg_area./cell_SA.*100;
difference = derived_cover - bergConc;

% Plot to check how we did
figure; 
subplot(1,3,1); imagesc(bergConc); colorbar; caxis([0 10]); title('target cover');
subplot(1,3,2); imagesc(derived_cover); colorbar; caxis([0 10]);  title('derived cover');
ax3= subplot(1,3,3); imagesc(difference); colorbar; caxis([-5 5]); title('difference');
print('-dpng','-r150','derived_iceberg_cover.png');
close

%% Now calculate iceberg surface area in every cell
% we will need to do this in MITgcm as well, but do it here as well to
% check and to plot

% pre-allocate
numBergsPerCell = zeros(size(bergMask)).*NaN;
total_long_face_SA = zeros([size(bergMask),length(depth)]).*NaN;
total_short_face_SA = zeros([size(bergMask),length(depth)]).*NaN;
total_base_SA = zeros([size(bergMask),length(depth)]).*NaN;
total_berg_SA = zeros([size(bergMask),length(depth)]).*NaN;
total_berg_volume = zeros([size(bergMask),length(depth)]).*NaN;

% Loop through cells and then through bergs within that cell, then through
% water column
disp('Calculating iceberg surface area in each cell...');
for i = 1:size(bergMask,1)
    for j = 1:size(bergMask,2)
        
        % Make sure we are in a cell containing icebergs
        if abs(bergMask(i,j)) == 1
            
            % Get lists of iceberg details
            local_berg_areas = all_berg_areas{i,j};
            local_berg_lengths = all_berg_lengths{i,j};
            local_berg_widths = all_berg_widths{i,j};
            local_berg_depths = all_berg_depths{i,j};
            
            % How many icebergs are there in this cell
            num_bergs = length(local_berg_depths);
            numBergsPerCell(i,j) = num_bergs; % keep track of number of bergs per cell (for MITgcm)
            
            % Check that we successfully put icebergs here
            if num_bergs == 0
                bergMask(i,j) = 0;
                numBergsPerCell(i,j) = 0;
                bergMaskNums(i,j) = 0;
            end
            
            % Pre-allocate arrays to hold temporary surface area
            % information (before computing totals)
            longface_SA = zeros([length(depth),length(num_bergs)]).*NaN;
            shortface_SA = zeros([length(depth),length(num_bergs)]).*NaN;
            base_SA = zeros([length(depth),length(num_bergs)]).*NaN;
            berg_volume = zeros([length(depth),length(num_bergs)]).*NaN;
            
            % Loop through icebergs
            for n = 1:num_bergs
                
                % get berg dimensions
                a = local_berg_areas(n);
                l = local_berg_lengths(n);
                w = local_berg_widths(n);
                d = local_berg_depths(n);
                
                % Get index of cell containing bottom of berg
                berg_bottom_cell = find(depth > d, 1, 'first');
                
                % Calculate surface area of individual faces and of all
                % faces in each cell
                if d < delz(1)
                    rem = d;
                else
                    rem = d - depth(berg_bottom_cell-1); % distance iceberg protrudes into bottom cell
                end
                longface_SA(1:berg_bottom_cell-1,n) = l.*squeeze(Z(i,j,1:berg_bottom_cell-1));
                longface_SA(berg_bottom_cell,n) = l.*rem;
                shortface_SA(1:berg_bottom_cell-1,n) = w.*squeeze(Z(i,j,1:berg_bottom_cell-1));
                shortface_SA(berg_bottom_cell,n) = w.*rem;
                base_SA(berg_bottom_cell,n) = a;
                
                % Calculate volume of each berg
                berg_volume(1:berg_bottom_cell-1,n) = l.*w.*squeeze(Z(i,j,1:berg_bottom_cell-1));
                berg_volume(berg_bottom_cell,n) = l.*w.*rem;
            end
                
            % Calculate total surface area at each depth
            total_long_face_SA(i,j,:) = nansum(shortface_SA.*2,2);
            total_short_face_SA(i,j,:) = nansum(longface_SA.*2,2);
            total_base_SA(i,j,:) = nansum(base_SA,2);
            total_berg_SA(i,j,:) =total_long_face_SA(i,j,:) + total_short_face_SA(i,j,:) + total_base_SA(i,j,:);
            total_berg_volume(i,j,:) = nansum(berg_volume,2);
            
        end
        
    end
end
 
% Calculate fraction of cell that is open
cell_volume = X.*Y.*Z; % m^3
cell_open_fraction = 1 - (total_berg_volume./cell_volume);

% Set NaNs to zero, for MITgcm
cell_open_fraction(isnan(cell_open_fraction)) = 0;
numBergsPerCell(isnan(numBergsPerCell)) = 0;


% Plot output?

%% Save
save('berg_generation_data.mat');
 
            
