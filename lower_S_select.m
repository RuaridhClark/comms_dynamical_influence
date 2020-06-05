function [trail] = lower_S_select(Scomp,n,in_list)
% Function that list all incoming node connections where the incoming node 
% has a smaller scalar projection.
    i = 1;
    trail(i) = n;
    for p = 1 : length(Scomp)
        if Scomp(p)>0
            i = i + 1;
            trail(i) = in_list{n}(p);
        end
    end
end