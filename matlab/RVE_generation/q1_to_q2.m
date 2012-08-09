%load coarse_mesh;
%p = [[0,1,2,0,0.7,1.3,2,0.7,1.3,0,1,2]',[0,0,0,1,0.7,0.7,1,1.3,1.3,2,2,2]'];
%t = [1,2,5; 2,6,5; 2,3,6; 3,7,6; 7,9,6; 7,12,9; 12,11,9; 11,8,9; 11,10,8; 10,4,8; 4,5,8; 4,1,5];
%nodes_free = logical([0,0,0,0,1,1,0,1,1,0,0,0]);
%[X,Y] = meshgrid(0:5,0:5); p = [X(:),Y(:)]; t = delaunay(p(:,1)',p(:,2)',{'Qt','Qz'});
%p = [0,0; 1,0; 1,1; 0,1; 0.5, 0.5]; t = [4,1,5; 1,2,5; 2,3,5; 3,4,5];
%nodes_free = logical(p(:,1)*0);

load([fname,'.mat'])

% Construct all edges
n = size(t,1);
edges = sort([t(:,1),t(:,2); t(:,2),t(:,3); t(:,3),t(:,1)],2);
et = [1:n,1:n,1:n]';
% Sorting makes it easier to detect duplicates
[tmp,I] = sort(edges(:,1)*max(edges(:,2)) + edges(:,2));
edges = edges(I,:);
et = et(I);

% Remove duplicates
et = [et,zeros(size(et))];
remove = zeros(size(edges,1),1);
for i = 2:size(edges,1)
    if edges(i,:) == edges(i-1,:)
        remove(i) = true;
        et(i-1,2) = et(i,1);
    end
end
%remove(I) = remove;
%edges(I,:) = edges;
%et(I,:) = et;
edges = edges(~remove,:);
et = et(~remove,:);

% Invert triangle->element list (et)
te = zeros(n,3);
k = zeros(n,1);
for i = 1:size(et,1)
    for j = 1:2
        q = et(i,j);
        if q > 0
            k(q) = k(q) + 1;
            te(q,k(q)) = i;
        end
    end
end

% Fix ordering... sigh, there must be some nicer way than all this shit.
for i = 1:size(t,1)
    ee = edges(te(i,:),:); % already sorted
    tt1 = sort(t(i,[1,2]));
    tt2 = sort(t(i,[2,3]));
    tt3 = sort(t(i,[3,1]));
    order = [find(ee(:,1) == tt1(1) & ee(:,2) == tt1(2)),...
             find(ee(:,1) == tt2(1) & ee(:,2) == tt2(2)),...
             find(ee(:,1) == tt3(1) & ee(:,2) == tt3(2))];
    te(i,:) = te(i,order);
end

% Add nodes
edges2 = [edges,(1:size(edges,1))'+size(p,1)];
new_p = (p(edges(:,1),:) + p(edges(:,2),:))/2;

% Construct new elements
p2 = [p;new_p];
t2 = [t,te+size(p,1)];

% Other stuff
nodes_left   = p(:,1) == min(p(:,1));
nodes_right  = p(:,1) == max(p(:,1));
nodes_top    = p(:,2) == max(p(:,2));
nodes_bottom = p(:,2) == min(p(:,2));

% For the corner nodes (needed to find boundaries only)
bc_nodes = nodes_left | nodes_right | nodes_top | nodes_bottom;

new_p_left   = new_p(:,1) < min(new_p(:,1))+1e-5;
new_p_right  = new_p(:,1) > max(new_p(:,1))-1e-5;
new_p_bottom = new_p(:,2) < min(new_p(:,2))+1e-5;
new_p_top    = new_p(:,2) > max(new_p(:,2))-1e-5;

nodes_left2   = p2(:,1) < min(p2(:,1))+1e-5;
nodes_right2  = p2(:,1) > max(p2(:,1))-1e-5;
nodes_bottom2 = p2(:,2) < min(p2(:,2))+1e-5;
nodes_top2    = p2(:,2) > max(p2(:,2))-1e-5;

bc_nodes2 = nodes_left2 | nodes_right2 | nodes_top2 | nodes_bottom2;

% Invert edges
ne = zeros(size(p,1),10);
k = zeros(size(p,1),1);
for i = 1:size(edges,1)
    for j = 1:2
        q = edges(i,j);
        k(q) = k(q) + 1;
        ne(q,k(q)) = i;
    end
end

% Find edges
x = ne(nodes_free,:); x = x(x>0);
edges_free = find(sparse(x,ones(size(x)),ones(size(x)))>1);

edges_left = find(new_p_left);
%x = ne(nodes_left,:); x = x(x>0);
%edges_left = find(sparse(x,ones(size(x)),ones(size(x)))>1);

edges_right = find(new_p_right);
%x = ne(nodes_right,:); x = x(x>0);
%edges_right = find(sparse(x,ones(size(x)),ones(size(x)))>1);

edges_top = find(new_p_top);
%x = ne(nodes_top,:); x = x(x>0);
%edges_top = find(sparse(x,ones(size(x)),ones(size(x)))>1);

edges_bottom = find(new_p_bottom);
%x = ne(nodes_bottom,:); x = x(x>0);
%edges_bottom = find(sparse(x,ones(size(x)),ones(size(x)))>1);
%

h = 0.01;
d = d_4particle_rve(new_p, 1.0, 1.05);
ix = edges_free; % New nodes have same numbering as edges.
nix = length(ix);
gradd = zeros(nix,2);
for ii=1:2
    a = zeros(1,2);
    a(ii) = h;
    d1x = d_4particle_rve(new_p(ix,:) + ones(nix,1)*a, 1.0, 1.05);
    gradd(:,ii)=(d1x-d(ix,:))/h;
end
new_p(ix,:) = new_p(ix,:) - d(ix)*ones(1,2).*gradd;

p2 = [p; new_p]; % Hack, I don't want to bother to rewrite this..

%return;
% Plot those bastards
hold on
%for i = 1:size(edges,1)
%    plot(p(edges(i,:),1),p(edges(i,:),2),'-k.')
%    plot(p2(edges2(i,[1,3,2]),1),p2(edges2(i,[1,3,2]),2),'-k.')
%end

%plot(new_p(:,1),new_p(:,2),'x')

plot(p2(nodes_left2,1),p2(nodes_left2,2),'bo')
plot(p2(nodes_right2,1),p2(nodes_right2,2),'bo')
plot(p2(nodes_top2,1),p2(nodes_top2,2),'ro')
plot(p2(nodes_bottom2,1),p2(nodes_bottom2,2),'ro')

for i = 1:length(edges_free)
    plot(p(edges(edges_free(i),:),1),p(edges(edges_free(i),:),2),'r');
end
for i = 1:length(edges_left)
    plot(p(edges(edges_left(i),:),1),p(edges(edges_left(i),:),2),'b','linewidth',2);
end
for i = 1:length(edges_right)
    plot(p(edges(edges_right(i),:),1),p(edges(edges_right(i),:),2),'b','linewidth',2);
end
for i = 1:length(edges_top)
    plot(p(edges(edges_top(i),:),1),p(edges(edges_top(i),:),2),'r','linewidth',2);
end
for i = 1:length(edges_bottom)
    plot(p(edges(edges_bottom(i),:),1),p(edges(edges_bottom(i),:),2),'r','linewidth',2);
end

%for i = 1:size(t2,1)
%    px = p2(t2(i,[1,4,2,5,3,6,1]),1);
%    py = p2(t2(i,[1,4,2,5,3,6,1]),2);
%    plot(px,py)
%end

axis equal
%axis ([-1.1,1.1,-1.1,1.1])

