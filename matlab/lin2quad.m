function [t2, p2] = lin2quad(t,p)
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

