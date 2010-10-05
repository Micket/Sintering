function m = v2m(v)
m = zeros(sqrt(length(v)));
if length(v) == 1
    order = [1];
elseif length(v) == 4
    order = [1,4,3,2];
else
    order = [1,5,9,4,8,7,2,6,3];
end
m(order) = v;
