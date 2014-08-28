load nodes_200_mat.txt
load pressure_200_mat.txt

nodes_pressure = [nodes_200_mat(:,1:2),pressure_200_mat];

save nodes_200.dat nodes_pressure -ASCII


