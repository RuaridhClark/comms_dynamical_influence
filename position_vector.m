function [S] = position_vector(x)
% Function that evaluates the magnitude of each node's position vector with
% respect to the origin of the Euclidean space.
    Nra = size(x,2);            % number of nodes
    S=zeros(Nra,1);
    for i = 1 : Nra
        S(i) = norm(x(:,i));    % Distance from origin
    end

end