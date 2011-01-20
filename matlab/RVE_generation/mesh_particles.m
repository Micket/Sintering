% Creates a mesh from particle data.

total = prod(res);

% Remove duplicate points
for x0 = 1:res(1)
for y0 = 1:res(2)
    if active(x0,y0)
        % Find neighbouring particles.
        foot_0 = [foot_x(x0,y0); foot_y(x0,y0)];
        b0 = foot_0 - 2*tubewidth;
        b1 = foot_0 + 2*tubewidth;
        ind0 = max(floor((b0-bb0)./dx),1);
        ind1 = min(ceil((b1-bb0)./dx),res');
        for x = ind0(1):ind1(1)
        for y = ind0(2):ind1(2)

            if (x ~= x0 || y ~= y0) && active(x,y)
                if norm(foot_0 - [foot_x(x,y);foot_y(x,y)]) < 0.01*dx(1);
                    active(x,y) = false;
                    if normal_x(x,y) == 0 && normal_y(x,y) == 0
                        normal_x(x0,y0) = 0;
                        normal_y(x0,y0) = 0;
                    end
                end
            end
        end
        end
    end
end
end

% Fetch all nodes;
non = sum(active(:));
n = 0;
nodes = zeros(non,length(res));
node_index = zeros(res);
for i = 1:total
    if active(i)
        n = n + 1;
        nodes(n,:) = [foot_x(i),foot_y(i)];
        node_index(i) = n;
    end
end

% Find edges
normal_projection = 0.5; % Limit for cutting off seperate paths.
noe = 0;
edges = zeros(5*non,length(res)); 
edge_length = zeros(5*non, 1);

for i = 1:total
    if active(i)
        % Find neighbouring particles.
        foot_0 = [foot_x(i); foot_y(i)];
        x0 = foot_0 - 2*tubewidth;
        x1 = foot_0 + 2*tubewidth;
        ind0 = max(floor((x0-bb0)./dx),1);
        ind1 = min(ceil((x1-bb0)./dx),res');
        n0 = [normal_x(i); normal_y(i)];
        
        c = 0;
        neighbour = [];
        dist = [];
        for x = ind0(1):ind1(1)
        for y = ind0(2):ind1(2)
            if active(x,y) && node_index(i) ~= node_index(x,y)
                n = [normal_x(x,y);normal_y(x,y)];
                if all(n0 == 0) || all(n == 0) || (n0'*n > normal_projection)
                    c = c + 1;
                    neighbour(c) = node_index(x,y);
                    dist(c) = norm(foot_0 - [foot_x(x,y);foot_y(x,y)]);
                end
            end
        end
        end

        % Sort these by shortest distance.
        [temp, I] = sort(dist);
        neighbour = neighbour(I);

        noe = noe + 1;
        edges(noe,:) = sort([node_index(i),neighbour(1)]);
        edges_length(noe) = dist(1);
        t1 = nodes(neighbour(1),:) - nodes(node_index(i),:);
        for j = 2:length(neighbour)
            t = nodes(neighbour(j),:) - nodes(node_index(i),:);
            if t*t1' < 0
                noe = noe + 1;
                edges(noe,:) = sort([node_index(i),neighbour(j)]);
                edges_length(noe) = dist(j);
                break;
            end
        end
    end
end

% Remove all duplicate edges.
edges = edges(1:noe,:);
% Sort the edges;
edges_id  = edges(:,1)  + edges(:,2)*1000;
[temp,I] = unique(edges_id);
edges = edges(I,:);
noe = size(edges,1);

clf
hold on
nodes_x = nodes(:,1);
nodes_y = nodes(:,2);
for i = 1:non
    text(nodes(i,1),nodes(i,2),num2str(i))
end
plot(nodes(:,1), nodes(:,2), 'rx')
plot(nodes_x(edges'), nodes_y(edges'),'b-');
for i = 1:noe
    text(mean(nodes_x(edges(i,:))),mean(nodes_y(edges(i,:))),['\color{red}',num2str(i)],'HorizontalAlignment','center')
end
hold off
axis equal
axis([bb0(1),bb1(1),bb0(2),bb1(2)]);
