function regions = regionize(vol, region_size, overlap, ...
                             fg_seed_vol, bg_seed_vol)
% break a volume into regions of a manageable size

    have_fg = 0;
    if nargin > 3
        have_fg = 1;
    end

    have_bg = 0;
    if nargin > 4
        have_bg = 1;
    end

    regions = {};

    rows = floor(1:region_size/overlap:size(vol, 1));
    for r1=rows
        r2 = min([r1 + region_size - 1, size(vol, 1)]);
        cols = floor(1:region_size/overlap:size(vol, 2));
        for c1=cols
            c2 = min([c1 + region_size - 1, size(vol, 2)]);
            slices = floor(1:region_size/overlap:size(vol, 3));
            for s1=slices
                s2 = min([s1 + region_size - 1, size(vol, 3)]);

                region.r1 = r1;
                region.r2 = r2;

                region.c1 = c1;
                region.c2 = c2;

                region.s1 = s1;
                region.s2 = s2;

                region.vol = vol(r1:r2, c1:c2, s1:s2);

                if have_fg
                    region.fg_seeds = fg_seed_vol(r1:r2, c1:c2, s1:s2);
                end

                if have_bg
                    region.bg_seeds = bg_seed_vol(r1:r2, c1:c2, s1:s2);
                end

                regions{end+1} = region;
            end
        end
    end
end
