function [in_list] = incoming_conns(M)
% Function that lists all incoming connections for each node.
    Nra = length(M);

    in_list = cell(Nra,1);
    for i = 1 : Nra
        for j = 1 : Nra
            if M(j,i) > 0
                in_list{i,1} = [in_list{i,1},j];
            end
        end
    end
end