load test.mat
load h20.mat

for i = 1:20
   Z(i) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1), 500000);
end
