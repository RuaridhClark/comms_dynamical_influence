function [x] = normalise_eigvec(V,incl)
% Function that normalises the eigenvectors and scales the first left
% eigenvector to have the joint largest entry when compared with the other
% input eigenvectors.
    n = length(incl);                           % number of input eigenvectors
    for i = 1 : n
        x(:,i) = real(V(:,incl(i)));            % only consider the real parts of eigenvectors
        x(:,i) = x(:,i)./sum(abs(x(:,i)));      % normalise eigenvector based on element magnitude
    end
    
    if n > 1 && incl(1)==1                      % Normalise first lef eigenvector relative to other eigenvectors
        mc = max(max(abs(x(:,2:end))));         % set first eigvector to have as large an element as other eigenvectors
        x(:,1) = x(:,1)*(mc/max(abs(x(:,1))));
    end
    x=x';
end