clear all; close all;
load test.mat

z=zeros(8,1);

testbatchdata = permute(testbatchdata, [1,3,2]);
testbatchdata = reshape(testbatchdata,size(testbatchdata,1) ...
    *size(testbatchdata,2),size(testbatchdata,3));

f = (sum(testbatchdata,1)+5)/size(testbatchdata,1);
parameter_bA_4_RTS = log(f) - log(1-f);

load h10.mat
z(1) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1), 100000)
% Z_10_SAMS = SAMS2(parameter_W, parameter_a, parameter_b, 1e-4)
% Z_10_TAP = TAP(parameter_W, parameter_b, parameter_a,1e-4)
% Z_10_RTS = RTS(parameter_W, parameter_a, parameter_bA_4_RTS, parameter_b, 0.25)
z(2) = sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(1)));

load h20.mat
z(3) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1), 400000)
% Z_20_SAMS = SAMS(parameter_W, parameter_b, parameter_a)
% Z_20_TAP = TAP(parameter_W, parameter_b, parameter_a)
% Z_20_RTS = RTS(parameter_W, parameter_a, parameter_b, 0.3)
z(4) = sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(3)));

load h100.mat
z(5) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1),500000)
% Z_100_SAMS = SAMS(parameter_W, parameter_b, parameter_a)
% Z_100_TAP = TAP(parameter_W, parameter_b, parameter_a)
% Z_100_RTS = RTS(parameter_W, parameter_a, parameter_b, 0.4)
z(6) = sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(5)));

load h500.mat
z(7) = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1),6000000)
% Z_500_SAMS = SAMS(parameter_W, parameter_b, parameter_a)
% Z_500_TAP = TAP(parameter_W, parameter_b, parameter_a)
% Z_500_RTS = RTS(parameter_W, parameter_a, parameter_b, 0.5)
z(8) = sum(exp(testbatchdata*parameter_b' + sum(log(1+exp(repmat(parameter_a,[size(testbatchdata,1),1])+ ...
    testbatchdata*parameter_W)),2)-z(7)));

save z z.mat