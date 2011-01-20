res = [40,40];
bb0 = [-1;-1]*1.2;
bb1 = [1;1]*1.2;
dx = (bb1-bb0)./(res'-1);
tubewidth = dx(1)*2;

active = logical(zeros(res));

min_dist = inf+zeros(res);

% Negative indicates inactive
foot_x = zeros(res);
foot_y = zeros(res);
%foot_z = -ones(res);

normal_x = zeros(res);
normal_y = zeros(res);
%normal_z = -ones(res);

[grid_x, grid_y] = meshgrid(linspace(bb0(1), bb1(1), res(1)), linspace(bb0(2), bb1(2), res(2)));
grid_x = grid_x';
grid_y = grid_y';

% Find the edge nodes (this only needs to be approximately correct, so i'll locate within a bounding box from the edge nodes;

% The edges of the RVE
lines = [-1,-1,-1,1;...
         -1,1,1,1;...
         1,1,1,-1;...
         1,-1,-1,-1];

for i = 1:size(lines,1)
    p0 = lines(i,1:2)';
    p1 = lines(i,3:4)';
    line_length = norm(p1-p0);
    v = (p1-p0)/line_length;
    % Find the bounding box for the segment
    x0 = min(a,b) - tubewidth;
    x1 = max(a,b) + tubewidth;
    n = [-v(2);v(1)];

    ind0 = max(floor((x0-bb0)./dx),1);
    ind1 = min(ceil((x1-bb0)./dx),res');

    for x = ind0(1):ind1(1)
    for y = ind0(2):ind1(2)
        g = [grid_x(x,y);grid_y(x,y)];
        % Now calculate the shortest distances;
        xi = (g-p0)'*v/line_length;
        edge_point = false;
        if xi < 0
            edge_point = true;
            xi = 0;
        elseif xi > 1
            edge_point = true;
            xi = 1;
        end
        % closest coordinate on line segment;
        c = a + v*line_length*xi;

        % Distance;
        dist = norm(g-c);

        if (dist < tubewidth)
            if (dist < min_dist(x,y))
                min_dist(x,y) = dist;
                foot_x(x,y) = c(1);
                foot_y(x,y) = c(2);
                if ~edge_point
                    normal_x(x,y) = n(1);
                    normal_y(x,y) = n(2);
                end
                active(x,y) = true;
            end
        end
    end
    end
end

% Define the radius
radius = 1.1;
angle = acos(1/radius);
circle_segment = [-1, 1, radius, 6/4*pi + angle, 2*pi - angle;...
                  1, 1, radius, pi + angle, 6/4*pi - angle;...
                  1, -1, radius, pi/2 + angle, pi - angle ;...
                  -1, -1, radius, angle, pi/2 - angle ];

for i = 1:size(circle_segment)
    a = circle_segment(i,1:2)';
    r = circle_segment(i,3);
    v0 = circle_segment(i,4);
    v1 = circle_segment(i,5);

    % Find the bounding box for the segment
    x0 = a - r - tubewidth;
    x1 = a + r + tubewidth;

    ind0 = max(floor((x0-bb0)./dx),1);
    ind1 = min(ceil((x1-bb0)./dx),res');

    for x = ind0(1):ind1(1)
    for y = ind0(2):ind1(2)
        g = [grid_x(x,y);grid_y(x,y)];
        
        w = atan2(g(2) - a(2), g(1) - a(1) );
        % Fix annoying piece of shit matlab
        w(w<0) = w + 2*pi;
        edge_point = false;
        if w < v0
            w = v0;
            edge_point = true;
        elseif w > v1
            w = v1;
            edge_point = true;
        end

        n = [cos(w);sin(w)];

        c = a + n*r;
        
        dist = norm(g-c);

        if (dist < tubewidth)
            if (dist < min_dist(x,y))
                min_dist(x,y) = dist;
                foot_x(x,y) = c(1);
                foot_y(x,y) = c(2);
                if ~edge_point  % Only for non-edgepoints
                    normal_x(x,y) = n(1);
                    normal_y(x,y) = n(2);
                end
                active(x,y) = true;
            end
        end
    end
    end
end
