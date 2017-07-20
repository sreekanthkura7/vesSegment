function [seg_vol, seg_prob] = segment_vessels_random_walker(vol, options)

% Segment an angiogram into foreground (vessels) and background. For
% efficiency, this method breaks the angiogram volume into overlapping
% sub-volumes and processes each in parallel, if desired.
%
% Parameters:
% - vol: 3D angiogram to segment
% - options: structure specifying options to the segmentation algorithm
%
% Option structure fields:
% - mask: 3D binary indicator of regions to process (not 0), and
%         regions to ignore (0)
%
% - region_size: [default 50] The number of voxels on a side for the
%                sub-volumes.
%
% - guess_seeds: Whether to automatically take a best guess at seed
%                locations based on analyzing image intensity. Either
%                this field, or the fg_seed_vol and bg_seed_vol should
%                be set.
%
% - fg_seed_vol: [default []] A volume the same size as vol, where
%                voxel values indicate whether a voxel is
%                identified as a foreground seed (1) or not (0)
%
% - bg_seed_vol: [default []] A volume the same size as vol, where
%                voxel values indicate whether a voxel is
%                identified as a background seed (1) or not (0)
%
% - fg_percentage: [default 1%] the percentage of the the highest
%                  valued voxels in each fg_win to assume are
%                  foreground seeds.
%
% - fg_win: [default 300] when guessing foreground seeds, break the
%           source volume into subvolumes with fg_win voxels on a
%           side.
%
% - bg_percentage: [default 50%] the percentage of the the lowest
%                  valued voxels in each bg_win to assume are
%                  background seeds.
%
% - bg_win: [default 75] when guessing background seeds, break the
%           source volume into subvolumes with bg_win voxels on a
%           side.
%
% - bounds: 3x2 matrix indicating the minimum and maximum angiogram
%           volume indices to be processed. Useful for recomputing an
%           segmentation only within a subset of the original volume.
%
% - verbose: [default 0] binary indicating whether extra info about
%            operation should be printed.
%
% Oliver Hinds <ohinds@nmr.mgh.harvard.edu> 2017-04-19
%
% TODO
% - guess the expected proportion of FG and BG seeds from
%   known vessel density and pixel size

    debug = 0;

    tic

    % argument handling
    if nargin < 1 || isempty(vol)
        error('vol is a required argument')
        return
    end

    if nargin < 2
        options = struct();
    end

    if ~isfield(options, 'fg_seeds') && ~isfield(options, 'guess_seeds')
        error('either fg_seeds or guess_seeds must be specified in options')
        return
    end

    if ~isfield(options, 'bg_seeds') && ~isfield(options, 'guess_seeds')
        error('either bg_seeds or guess_seeds must be specified in options')
        return
    end

    if ~isfield(options, 'mask') || isempty(options.mask)
        options.mask = ones(size(vol));
    end
    mask = options.mask;

    if ~isfield(options, 'guess_seeds')
        options.guess_seeds = 0;
    end
    guess_seeds = options.guess_seeds;

    if ~isfield(options, 'region_size')
        options.region_size = 50;
    end
    region_size = options.region_size;

    if ~isfield(options, 'fg_seed_vol')
        options.fg_seed_vol = [];
    end
    fg_seed_vol = options.fg_seed_vol;

    if ~isfield(options, 'bg_seed_vol')
        options.bg_seed_vol = [];
    end
    bg_seed_vol = options.bg_seed_vol;

    if ~isfield(options, 'fg_percentage')
        options.fg_percentage = 1;
    end
    fg_percentage = options.fg_percentage;

    if ~isfield(options, 'fg_win')
        options.fg_win = 300;
    end
    fg_win = options.fg_win;

    if ~isfield(options, 'bg_percentage')
        options.bg_percentage = 50;
    end
    bg_percentage = options.bg_percentage;

    if ~isfield(options, 'bg_win')
        options.bg_win = 50;
    end
    bg_win = options.bg_win;

    if ~isfield(options, 'bounds') || isempty(options.bounds)
        whole_vol = 1;
    else
        whole_vol = 0;
        ret_seg_vol = zeros(size(vol));
        ret_seg_prob = zeros(size(vol));
        vol = vol(options.bounds(1,1):options.bounds(1,2), ...
                  options.bounds(2,1):options.bounds(2,2), ...
                  options.bounds(3,1):options.bounds(3,2));
        fg_seed_vol = fg_seed_vol(options.bounds(1,1):options.bounds(1,2), ...
                                  options.bounds(2,1):options.bounds(2,2), ...
                                  options.bounds(3,1):options.bounds(3,2));
        bg_seed_vol = bg_seed_vol(options.bounds(1,1):options.bounds(1,2), ...
                                  options.bounds(2,1):options.bounds(2,2), ...
                                  options.bounds(3,1):options.bounds(3,2));
        bounds = options.bounds;
    end

    if ~isfield(options, 'verbose')
        options.verbose = 0;
    end
    verbose = options.verbose;

    clear options;

    % prepare outputs
    seg_vol = zeros(size(vol));
    seg_prob = zeros(size(vol));

    % find the seeds
    if guess_seeds
        fg_seed_vol = find_seeds_by_threshold(...
            vol, fg_percentage, 1, fg_win, [], verbose);
        bg_seed_vol = find_seeds_by_threshold(...
            vol, bg_percentage, 0, bg_win, [], verbose);
    end

    fg_seed_vol(find(mask == 0)) = 0;
    bg_seed_vol(find(mask == 0)) = 1;

    % split the volume into parallelizable regions
    regions = regionize(vol, region_size, 1.5, fg_seed_vol, bg_seed_vol);

    % segment each region
    num_regions = length(regions);
    parfor region_ind=1:length(regions)
        if verbose
            fprintf('%d/%d\n', region_ind, num_regions);
        end
        regions{region_ind} = segment_region(regions{region_ind}, ...
                                             verbose, debug);
    end

    % recombine the regions into a single volume
    for region_ind=1:length(regions)
        region = regions{region_ind};

        seg_vol(region.r1:region.r2, ...
                region.c1:region.c2, ...
                region.s1:region.s2) = ...
            seg_vol(region.r1:region.r2, ...
                    region.c1:region.c2, ...
                    region.s1:region.s2) | region.seg_region;

        seg_prob(region.r1:region.r2, ...
                 region.c1:region.c2, ...
                 region.s1:region.s2) = ...
            max(seg_prob(region.r1:region.r2, ...
                         region.c1:region.c2, ...
                         region.s1:region.s2), region.seg_prob);
    end

    if debug
        show_vol_seg(vol, seg_vol);
    end

    % if we only did part of the volume, put that part in the right place
    if ~whole_vol
        ret_seg_vol(bounds(1,1):bounds(1,2),...
                    bounds(2,1):bounds(2,2),...
                    bounds(3,1):bounds(3,2)) = seg_vol;
        seg_vol = ret_seg_vol;

        ret_seg_prob(bounds(1,1):bounds(1,2),...
                     bounds(2,1):bounds(2,2),...
                     bounds(3,1):bounds(3,2)) = seg_prob;
        seg_prob = ret_seg_prob;
    end

    toc
end

function region = segment_region(region, verbose, debug)

    if sum(region.vol(:) ~= 0) < 0.1 * prod(size(region.vol))

        if verbose
            fprintf(...
                'skipping mostly blank (%d / %d) region at %d %d %d\n', ...
                sum(region.vol(:) ~= 0), prod(size(region.vol)), ...
                region.r1, region.c1, region.s1);
        end

        region.seg_region = zeros(size(region.vol));
        region.seg_prob = zeros(size(region.vol));

        return;
    end

    fg_seeds = find(region.fg_seeds);
    bg_seeds = find(region.bg_seeds);

    if length(fg_seeds) < size(region.vol, 1)
        region.seg_region = zeros(size(region.vol));
        region.seg_prob = zeros(size(region.vol));
        return;
    end

    if length(bg_seeds) < size(region.vol, 1)
        region.seg_region = zeros(size(region.vol));
        region.seg_prob = zeros(size(region.vol));
        return;
    end

    if verbose
        fprintf(['segmenting a volume at %d,%d,%d of size %dx%dx%d ' ...
                 'with %d fg and %d bg seeds\n'], region.r1, region.c1, ...
                region.s1, size(region.vol), length(fg_seeds), ...
                length(bg_seeds));
    end

    [seg_region seg_prob] = random_walker_3d(...
        region.vol, [fg_seeds; bg_seeds], ...
        [ones(1, size(fg_seeds, 1)), ...
         2 * ones(1, size(bg_seeds, 1))]);
    seg_region = seg_region == 1;
    seg_prob = seg_prob(:,:,:,1);

    region.seg_region = seg_region;
    region.seg_prob = seg_prob;

    if debug
        h = show_vol_seg(region.vol, region.seg_region);
        pause;
        close all;
    end
end
