function [newVals] = phaseShift(vals, shift, loc_nada) % takes shift in degrees
    shift = mod(shift, 360); % mods phase shift by 360
    newVals = zeros(size(vals)); %New array of values pre-allocated
    period = loc_nada(3) - loc_nada(1); % period determined by zero times
    cut_off = round(period*(shift/360));
    newVals(1:end-cut_off) = vals(cut_off+1:end);
    newVals(end-cut_off:end) = -vals(end:-1:end-cut_off);
    
end