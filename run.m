clear all;
load h10.mat

Z_10_TAP = TAP(parameter_W, parameter_b, parameter_a)
Z_10_AIS = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1))

% load test.mat
% testbatchdata = reshape(permute(testbatchdata, [1,3,2]),size(testbatchdata,1) ...
%     *size(testbatchdata,2),size(testbatchdata,3));
% 