clear all; close all;
load test.mat

z=zeros(8,1);

testbatchdata = permute(testbatchdata, [1,3,2]);
testbatchdata = reshape(testbatchdata,size(testbatchdata,1) ...
    *size(testbatchdata,2),size(testbatchdata,3));

load h10.mat
z(1) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1), 30000);
z(5) = log(sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(1))));

load h20.mat
z(2) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1), 200000);
z(6) = log(sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(2))));

load h100.mat
z(3) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1),60000);
z(7)= log(sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(3))));

load h500.mat
z(4) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1),60000);
z(8) = log(sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(4))));

save z.mat z