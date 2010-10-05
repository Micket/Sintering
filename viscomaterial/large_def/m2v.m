function v = m2v(m)
if size(m,1) == 1
    order = [1];
elseif size(m,1) == 2
    order = [1,4,3,2];
else
    order = [1,5,9,4,8,7,2,6,3];
end
v = m(order)';

