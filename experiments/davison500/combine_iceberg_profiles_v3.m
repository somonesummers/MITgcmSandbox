% combine_iceberg_profiles.m
% script to support iceberg package in MITgcm
% It searches through MITgcm run directory for profiles of iceberg melt
% rate, freshwater flux, heat flux and current speed relative to iceberg
% drift velocity. It combines these into a 3D matrix (z, berg, time), with
%  corresponding vectors of berg dimensions and location within
% the domain (x,y)

clear
close all
clc

% Load berg mask and mask nums
load('berg_generation_data.mat','bergMaskNums','nz','numBergsPerCell');

% Get list of text files (e.g. iceberg_data_00001_000001_00000001.txt
% i.e. iceberg_data_cellNum_bergNum_myIter_counter.txt)
fileList_melt = dir('iceberg_mltRt*.txt');
fileList_fw = dir('iceberg_fwFlx*.txt');
fileList_hf = dir('iceberg_htFlx*.txt');
fileList_sp1 = dir('iceberg_spdla*.txt');
fileList_sp2 = dir('iceberg_spdsa*.txt');
fileList_sp3 = dir('iceberg_spdba*.txt');

% pre allocate big arrays to hold output (trim later)
mProf = zeros([nz,nansum(numBergsPerCell(:)),100]).*NaN; % max 100 dumps?
fwProf = zeros([nz,nansum(numBergsPerCell(:)),100]).*NaN; % max 100 dumps?
hfProf = zeros([nz,nansum(numBergsPerCell(:)),100]).*NaN; % max 100 dumps?
rSp1Prof = zeros([nz,nansum(numBergsPerCell(:)),100]).*NaN; % max 100 dumps?
rSp2Prof = zeros([nz,nansum(numBergsPerCell(:)),100]).*NaN; % max 100 dumps?
rSp3Prof = zeros([nz,nansum(numBergsPerCell(:)),100]).*NaN; % max 100 dumps?
x = zeros([nansum(numBergsPerCell(:)),1]).*NaN;
y = zeros([nansum(numBergsPerCell(:)),1]).*NaN;
depth_data = zeros([nansum(numBergsPerCell(:)),1]).*NaN;
length_data = zeros([nansum(numBergsPerCell(:)),1]).*NaN;
width_data = zeros([nansum(numBergsPerCell(:)),1]).*NaN;


% Read data
count = 0;
for i = 1:numel(fileList_melt)
    % load data
    fid = fopen(fileList_melt(i).name);
    melt_data = textscan(fid,'%s%s%s%s','CollectOutput',true); 
    fclose(fid);
    fid = fopen(fileList_fw(i).name); fw_data = textscan(fid,'%s%s%s%s','CollectOutput',true); fclose(fid);
    fid = fopen(fileList_hf(i).name); hf_data = textscan(fid,'%s%s%s%s','CollectOutput',true); fclose(fid);
    fid = fopen(fileList_sp1(i).name); sp1_data = textscan(fid,'%s%s%s%s','CollectOutput',true); fclose(fid);
    fid = fopen(fileList_sp2(i).name); sp2_data = textscan(fid,'%s%s%s%s','CollectOutput',true); fclose(fid);
    fid = fopen(fileList_sp3(i).name); sp3_data = textscan(fid,'%s%s%s%s','CollectOutput',true); fclose(fid);
    melt_data = melt_data{:}';
    melt_data = melt_data(:);
    fw_data = fw_data{:}'; fw_data = fw_data(:);
    hf_data = hf_data{:}'; hf_data = hf_data(:);
    sp1_data = sp1_data{:}'; sp1_data = sp1_data(:);
    sp2_data = sp2_data{:}'; sp2_data = sp2_data(:);
    sp3_data = sp3_data{:}'; sp3_data = sp3_data(:);
    % extract data
    cellNum = str2double(cell2mat(melt_data(6)));
    [row,col] = find(bergMaskNums == cellNum);
    numBergs = numBergsPerCell(row,col);
    melt_rate_cell = melt_data(1:8:end);
    berg_depth_cell = melt_data(2:8:end);
    berg_length_cell = melt_data(3:8:end);
    berg_width_cell = melt_data(5:8:end);
    iter_cell = melt_data(7:8:end);
    zlevel_cell = melt_data(8:8:end);
    fw_cell = fw_data(1:8:end);
    hf_cell = hf_data(1:8:end);
    sp1_cell = sp1_data(1:8:end);
    sp2_cell = sp2_data(1:8:end);
    sp3_cell = sp3_data(1:8:end);
    % pre-allocate double arrays
    melt_rate_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    berg_depth_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    berg_length_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    berg_width_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    iter_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    zlevel_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    fw_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    hf_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    sp1_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    sp2_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    sp3_dbl = zeros([numel(melt_rate_cell),1]).*NaN;
    % convert to double
    for j = 1:numel(melt_rate_cell)
        dummy = str2double(cell2mat(melt_rate_cell(j)));
        melt_rate_dbl(j) = dummy;
        dummy = str2double(cell2mat(berg_depth_cell(j)));
        berg_depth_dbl(j) = dummy;
        dummy = str2double(cell2mat(berg_length_cell(j)));
        berg_length_dbl(j) = dummy;
        dummy = str2double(cell2mat(berg_width_cell(j)));
        berg_width_dbl(j) = dummy;
        dummy = str2double(cell2mat(iter_cell(j)));
        iter_dbl(j) = dummy;
        dummy = str2double(cell2mat(zlevel_cell(j)));
        zlevel_dbl(j) = dummy;
        dummy = str2double(cell2mat(fw_cell(j)));
        fw_dbl(j) = dummy;
        dummy = str2double(cell2mat(hf_cell(j)));
        hf_dbl(j) = dummy;
        dummy = str2double(cell2mat(sp1_cell(j)));
        sp1_dbl(j) = dummy;
        dummy = str2double(cell2mat(sp2_cell(j)));
        sp2_dbl(j) = dummy;
        dummy = str2double(cell2mat(sp3_cell(j)));
        sp3_dbl(j) = dummy;
        clear dummy
    end
    % Iterations
    uniqueIter = unique(iter_dbl(:));
    numIters = length(uniqueIter);
    % convert to desired shape
    melt_data_tmp = reshape(melt_rate_dbl,[nz,numBergs,numIters]);
    depth_data_tmp = reshape(berg_depth_dbl,[nz,numBergs,numIters]);
    length_data_tmp = reshape(berg_length_dbl,[nz,numBergs,numIters]);
    width_data_tmp = reshape(berg_width_dbl,[nz,numBergs,numIters]);
    fw_data_tmp = reshape(fw_dbl,[nz,numBergs,numIters]);
    hf_data_tmp = reshape(hf_dbl,[nz,numBergs,numIters]);
    sp1_data_tmp = reshape(sp1_dbl,[nz,numBergs,numIters]);
    sp2_data_tmp = reshape(sp2_dbl,[nz,numBergs,numIters]);
    sp3_data_tmp = reshape(sp3_dbl,[nz,numBergs,numIters]);
    % prepare bergs index
    start_idx = count+1;
    end_idx = count+numBergs;
    % Send to large arrays
    mProf(:,start_idx:end_idx,1:numIters) = melt_data_tmp;
    fwProf(:,start_idx:end_idx,1:numIters) = fw_data_tmp;
    hfProf(:,start_idx:end_idx,1:numIters) = hf_data_tmp;
    rSp1Prof(:,start_idx:end_idx,1:numIters) = sp1_data_tmp;
    rSp2Prof(:,start_idx:end_idx,1:numIters) = sp2_data_tmp;
    rSp3Prof(:,start_idx:end_idx,1:numIters) = sp3_data_tmp;
    % Get dimensions and locations of bergs
    bergDepth(start_idx:end_idx) = squeeze(depth_data_tmp(1,:,1));
    bergWidth(start_idx:end_idx) = squeeze(width_data_tmp(1,:,1));
    bergLength(start_idx:end_idx) = squeeze(length_data_tmp(1,:,1));
    x(start_idx:end_idx) = repmat(row,[length(start_idx:end_idx),1]);
    y(start_idx:end_idx) = repmat(col,[length(start_idx:end_idx),1]);
    % Increment count (how many icebergs done so far)
    count = end_idx;
    % clear variables
    clear fid melt_data fw_data hf_data sp1_data sp2_data sp3_data cellNum row col numBergs
    clear melt_rate_cell berg_depth_cell berg_width_cell berg_length_cell iter_cell zlevel_cell
    clear fw_cell hf_cell sp1_cell sp2_cell sp3_cell
    clear uniqueIter numIters melt_data_tmp depth_data_tmp length_data_tmp width_data_tmp
    clear fw_data_tmp hf_data_tmp sp1_data_tmp sp2_data_tmp sp3_data_tmp start_idx end_idx
end

    
save('iceberg_profile_data.mat','mProf','fwProf','hfProf','rSp1Prof','rSp2Prof','rSp3Prof','bergDepth','bergLength','bergWidth','x','y');
