% N = 50;  % Number of points
% x_pos = rand(N, 1);  % X-coordinates between 0 and 1
% y_pos = rand(N, 1);  % Y-coordinates between 0 and 1
%  
% distance_matrix = pdist2([x_pos, y_pos], [x_pos, y_pos]);
% 
% k = 7;  % Number of nearest neighbors
% [idx, D] = knnsearch([x_pos, y_pos], [x_pos, y_pos], 'K', k+1);  % +1 to include the point itself
% 
% adj = zeros(N, N);
% for i = 1:N   % fix this
%     adj(i, idx(i, 1:end)) = 1;  % Exclude the point itself
% end
% adj = adj - diag(diag(adj));

% [C,x,S_vals] = CDI_nopath(adj,'A',[1,2,3]);

L = diag(diag(adj))-adj;

[C,x,S_vals] = CDI_nopath(L,'L',[1,2,3]);

DG =digraph(adj);

figure
for i = 1 : max(C)
    Inds = find(C==i);
    scatter3(x(Inds,1),x(Inds,2),x(Inds,3),0.01+(abs(S_vals(Inds,:)))*1000,'filled')
    hold on
end
xlabel('v_1','fontsize', 14)
ylabel('v_2','fontsize', 14)
zlabel('v_3','fontsize', 14)

%%
figure
LWidths = 1*DG.Edges.Weight/max(DG.Edges.Weight);
p = plot(DG,'LineWidth',LWidths,'EdgeColor',[.5,.5,.5],'NodeColor',[.5,.5,.5],'EdgeAlpha',.15);
% To remove all labels
p.NodeLabel = {};
p.XData = x(:,1);
p.YData = x(:,2);
p.ZData = x(:,3);
hold on
scatter3(x(:,1),x(:,2),x(:,3),0.01+(abs(S_vals))*1000,'filled','LineWidth',2)

%%
figure
LWidths = 1*DG.Edges.Weight/max(DG.Edges.Weight);
p = plot(DG,'LineWidth',LWidths,'EdgeColor',[.5,.5,.5],'NodeColor',[.5,.5,.5]);
% To remove all labels
p.NodeLabel = {};
p.XData = x_pos;
p.YData = y_pos;
axis off

%%
figure
LWidths = 1*DG.Edges.Weight/max(DG.Edges.Weight);
p = plot(DG,'LineWidth',LWidths,'EdgeColor',[.5,.5,.5],'NodeColor',[.5,.5,.5],'EdgeAlpha',.5);
% To remove all labels
p.NodeLabel = {};
p.XData = x_pos;
p.YData = y_pos;
hold on
scatter(x_pos,y_pos,0.01+(abs(S_vals))*1000,'filled','LineWidth',2)
axis off