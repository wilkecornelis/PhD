clear all
clc

l_range = linspace(-1,1,1000);
lambda = 0.15;
k = 2*pi/lambda;

d = 0.15;
N_elements = 10;

d_vec = linspace(1,N_elements*0.150,N_elements);
g_vec = ones(N_elements,1);

g_vec(5) = 1;
g_vec(8) = 0;
g_vec(6) = 1;

% s_points = linspace(500-(N_elements/2)+1,500+(N_elements/2),N_elements);
s_points = linspace(1,1000,N_elements);
% s_points = linspace(420,580,N_elements);
% 

for n = 1:length(l_range)
    AF(n) = sum(g_vec.*transpose(exp(j*k.*d_vec*l_range(n))));
    
end    

b_vec = AF(round(s_points));

for n = 1:N_elements
    A_mat(:,n) = transpose(exp(j*k.*d_vec*l_range(round(s_points(n)))));

end

A_mat_inv = inv(A_mat);
b_vec_tr = transpose(b_vec);

b_vec_solve = b_vec*A_mat_inv



