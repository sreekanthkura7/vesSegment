function seed_vol = find_seeds_by_threshold(...
    vol, percentage, fg, region_size, bounds, verbose)

    debug = 0;

    if verbose
        tic
    end

    if nargin < 5 || isempty(bounds)
        whole_vol = 1;
    else
        whole_vol = 0;
        orig_vol_size = size(vol);
        vol = vol(bounds(1,1):bounds(1,2), ...
                  bounds(2,1):bounds(2,2), ...
                  bounds(3,1):bounds(3,2));
    end

    regions = regionize(vol, region_size, 1);

    num_regions = length(regions);
    parfor region_ind=1:length(regions)
        if verbose
            fprintf('%d/%d\n', region_ind, num_regions);
        end
        regions{region_ind} = find_seeds_region(...
            regions{region_ind}, percentage, fg, verbose, debug);
    end

    seed_vol = zeros(size(vol));
    for region_ind=1:length(regions)
        region = regions{region_ind};

        seed_vol(region.r1:region.r2, ...
                 region.c1:region.c2, ...
                 region.s1:region.s2) = ...
                   region.seed_vol;
    end

    if ~whole_vol
        ret_seed_vol = zeros(orig_vol_size);

        ret_seg_vol(bounds(1,1):bounds(1,2),...
                    bounds(2,1):bounds(2,2),...
                    bounds(3,1):bounds(3,2)) = seed_vol;
        seed_vol = ret_seg_vol;
    end

    if verbose
        toc
    end
end

function region = find_seeds_region(region, percentage, fg, verbose, debug)

% guess seeds as % of the brightest (or darkest) pixels
    if fg
        sorted_vox_vals = sort(region.vol(:), 'descend');
        thresh = sorted_vox_vals(...
            round(percentage / 100 * length(sorted_vox_vals)));
        seeds = find(region.vol(:) > thresh);
    else
        sorted_vox_vals = sort(region.vol(:), 'ascend');
        thresh = sorted_vox_vals(...
            round(percentage / 100 * length(sorted_vox_vals)));
        seeds = find(region.vol(:) <= thresh);
    end

    region.seed_vol = zeros(size(region.vol));
    region.seed_vol(seeds) = 1;
end
