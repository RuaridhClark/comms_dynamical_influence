function [C,ev_vals,S_vals] = CDI(M,Type,Incl)
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

    M_sparse=sparse(M);             % sparse matrix
    Nra= length(M);                 % length of Adj
    if Type == 'L'                  % if using Laplacian matrix
        M=-M;           
    end
    [V,D]=eig(M');                  % calculates the left eigenvectors
    [~,J]=sort(real(diag(D)),'desc');
    V=V(:,J);                       % order according to eigenvalues from largest to smallest
    [x] = normalise_eigvec(V,Incl); % normalise and scale eigenvectors
    [S] = position_vector(x);     % calculate position vector for each node in eigenvector space
    [~,I] = sort(S,'descend');      % sort nodes according to vector magnitude
    [in_list] = incoming_conns(M);  % Create a list of incoming connections for each node

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
    
    group=cell(1,length(leaders)); assigned=[]; 
    for k = 1 : length(leaders)
        cn = leaders(k); 
        rmvd_list = [];
        node_list = cn;             % create a list of nodes whose connections have not been searched
        S_cn = S(cn);               % magnitude of Euclidean position vector for node cn
        if ~isempty(in_list{cn})    % if node cn has any incoming connections then continue
            Sproj=Ssave{k};         
            while isempty(node_list)==0 % while there are nodes still to search continue
                fn=node_list(1);        % next node on node_list
                if ~ismember(fn,assigned)
                    viable=k : length(leaders); % fn nodes can only be assigned to communities not yet defined
                    [~,mI]=sort(dot(x(:,leaders(viable)),repmat(x(:,fn),1,length(leaders(viable))))./(S(leaders(viable))'),'desc');
                    if leaders(viable(mI(1)))~=leaders(k)
                        for i = 1 : length(mI)
                            [~,path,~] = graphshortestpath(M_sparse,fn,leaders(viable(mI(i))),'Method','BFS'); % MATLAB in-built function for shortest path using Breadth first search (BFS)
                            if ~isempty(path)
                                mI=mI(i);       % selected leader with closest alignment AND a path from fn to leader
                                break;
                            end
                        end
                    end
                    if leaders(viable(mI(1))) == leaders(k)
                        Vt = x(:,in_list{fn});                                  % position vectors for all fn's incoming nodes
                        [Sproj] = scalar_proj(x(:,cn),Vt,S_cn);                 % dot product normalised by magnitude cn's position vector 
                        [Sdot_fn] = scalar_proj(x(:,cn),x(:,fn),S_cn);          % fn position vector magnitude in cn direction
                        Scomp = Sdot_fn - Sproj;                                % Compare with the in_list nodes
                        [trail] = lower_S_select(Scomp,fn,in_list);             % select all incoming nodes to fn that have a smaller Sdot than fn
                        trail_red=trail(~ismember(trail,group{k}));             % trail does not include nodes already assigned
                        group{k}=[group{k},fn];
                        assigned=[assigned,fn];                                 % list of nodes assigned to a community
                        sdot_g{k}=[sdot_g{k},Sproj(ismember(in_list{fn},trail_red))]; % all sdot values for each group wrt to cn
                        trail_red=trail_red(~ismember(trail_red,node_list));    % trail does not include nodes on node_list
                        trail_red=trail_red(~ismember(trail_red,rmvd_list));    % trail does not include node already searched for this community
                        node_list=[node_list,trail_red];
                    end
                end
                rmvd_list = [rmvd_list,node_list(1)];
                node_list(1)=[];
            end
        end
    end

    ev_vals=zeros(length(M),length(group));
    S_vals=zeros(length(M),length(group));
    for i = 1 : length(group)
        g = group{i};
        for j= 1 : length(g)
            ev_vals(g(j),i)=x(1,g(j)); 	% first eigenvector values
            S_vals(g(j),i)=S(g(j));     % position vector values
        end
    end

    mv = max(abs(ev_vals));             % find max first eigenvector in each community
    [~,Imv] = sort(mv,'desc');          % sort according to max first eigenvector
    
    C = zeros(Nra,1); 
    vecs_ind=abs(ev_vals)>0;
    for j = 1 : length(Imv)
        C(vecs_ind(:,Imv(j)))=j;        % Community designation where order is determined by first eigenvector          
    end 
    ev_vals=ev_vals(:,Imv);             % reordered to match community order
end