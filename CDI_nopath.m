function [C,x,S_vals] = CDI_nopath(M,Type,Incl)
% CDI divides a network into communities of dynamical influence.
%
% [C,ev_vals,S_vals] = CDI(M,type,incl) determines the communities of 
% of dynamical influence for a matrix M of either Type:
% 'A' - adjacency matrix
% 'L' - Laplacian matrix
% The number of eigenvectors included in the evaluation (typically 3) is
% defined by Incl ([1,2,3] would include the first 3 eigenvectors of matrix
% M in CDI evaluation.
% Outputs:  C       - vector listing each node's community designation
%           ev_vals - matrix of first eigenvalue entries where each column 
%                     contains contains nodes from a single community
%           S_vals  - matrix of first eigenvalue entries where each column 
%                     contains contains nodes from a single community
%
% References:   Clark, R., Punzo, G. & Macdonald, M. Network Communities of
%               Dynamical Influence. Sci Rep 9, 17590 (2019).

    [G,Nra,M] = graph_setup(M,Type);
    
    [S,x,I,in_list] = eig_embedding(M,Incl);

    [leaders,sdot_g,Ssave] = find_leaders(S,x,Nra,I,in_list); 

    C = zeros(Nra,1); 
    S_vals = zeros(Nra,length(leaders));
    for i = 1: length(S)
        [Sproj] = scalar_proj(x(:,leaders),x(:,i),S(i));                 % dot product normalised by magnitude cn's position vector 
        [~,I_max] = max(Sproj);
        C(i)=I_max;
        S_vals(i,C(i)) = S(i);
    end
    x=x';
end

%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%
function [G,Nra,M] = graph_setup(M,Type)
    G = digraph(M);
    Nra= length(M);                 % length of Adj
    if Type == 'L'                  % if using Laplacian matrix
        M=-M;           
    end
end

function [S,x,I,in_list] = eig_embedding(M,Incl)
%     [V,D]=eig(M');                  % calculates the left eigenvectors

    num_exclude = Incl(1); % Number of largest eigenvalues/eigenvectors to exclude
    num_eig = Incl(end); % Total number of eigenvalues/eigenvectors to compute
    
    opts = struct('largest', num_eig, 'smallest', num_exclude);
    [V, D] = eigs(M', num_eig, 'lm', opts);
    
    V = V(:, num_exclude:end); % Select eigenvectors from num_exclude+1 to the end
    D = D(:, num_exclude:end); % Select eigenvectors from num_exclude+1 to the end

    [~,J]=sort(real(diag(D)),'desc');
    V=V(:,J);                       % order according to eigenvalues from largest to smallest
    [x] = normalise_eigvec(V,1:length(Incl)); % normalise and scale eigenvectors
    [S] = position_vector(x);     % calculate position vector for each node in eigenvector space
    [~,I] = sort(S,'descend');      % sort nodes according to vector magnitude
    [in_list] = incoming_conns(M);  % Create a list of incoming connections for each node
end

function [leaders,sdot_g,Ssave] = find_leaders(S,x,Nra,I,in_list)
    k=0; group={};
    for i = 1 : Nra
        cn = I(i);                  % Identify the node with the i-th largest euclidean position vector 
        S_cn = S(cn);               % Magnitude of Euclidean position vector for node cn

        if ~isempty(in_list{cn})    % sum(M(:,cn))>0 % if node cn has any incoming connections then continue
            Vt = x(:,in_list{cn});  % position vector for all incoming connections
            [Sproj] = scalar_proj(x(:,cn),Vt,S_cn); %Dot product for determining a nodes S value with respect to the PN
            if (min(S_cn-Sproj)>=0)  % Check that all of cn's incoming connections have either smaller FLE or S value
                k=k+1;
                leaders(k)=cn;      % List of each community's nodes
                sdot_g{k}=S_cn;     % List of Sdots wrt cn
                Ssave(k)={Sproj};   % Save Sdot for later use
            end
        end
    end
end

function [mI] = path2leader(G,fn,mI,leaders,viable,k)
    if leaders(viable(mI(1)))~=leaders(k)
        for i = 1 : length(mI)
            [path,~,~] = shortestpath(G,fn,leaders(viable(mI(i))),'Method','auto');
%                             [~,path,~] = graphshortestpath(M_sparse,fn,leaders(viable(mI(i))),'Method','BFS'); % MATLAB in-built function for shortest path using Breadth first search (BFS)
            if ~isempty(path)
                mI=mI(i);       % selected leader with closest alignment AND a path from fn to leader
                break;
            end
        end
    end
end

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

function [Sdot] = scalar_proj(vec_a,vec_b,magn)

    if size(vec_a,2)==1 && size(vec_b,2)~=1
        vec_a=repmat(vec_a,1,size(vec_b,2));
    elseif size(vec_b,2)==1 && size(vec_a,2)~=1
        vec_b=repmat(vec_b,1,size(vec_a,2));
    end

    Sdot = dot(vec_a,vec_b)./(magn);

end