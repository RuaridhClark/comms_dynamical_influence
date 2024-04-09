function [S] = position_vector(x)
% Function that evaluates the magnitude of each node's position vector with
% respect to the origin of the Euclidean space.
    S=vecnorm(x,2);

end