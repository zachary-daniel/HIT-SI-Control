
function [ts] = locsToTimes(locs, time)
    for i = 1:length(locs)
        t = locs(i);
        locs(i) = time(t);
    end
    ts = locs;
end
