clf
hold on

for i = 1:numel(active)
if active(i)
    plot([grid_x(i),foot_x(i)],[grid_y(i),foot_y(i)],'c');
    plot(grid_x(i), grid_y(i), 'go');
    plot(foot_x(i), foot_y(i), 'rx');
    quiver(foot_x(i), foot_y(i), normal_x(i), normal_y(i), 0.1);
end
end
axis equal
