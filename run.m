clear all; close all;
load test.mat
testbatchdata = permute(testbatchdata, [1,3,2]);
testbatchdata = reshape(testbatchdata,size(testbatchdata,1) ...
    *size(testbatchdata,2),size(testbatchdata,3));

f = (sum(testbatchdata,1)+5)/size(testbatchdata,1);
parameter_bA_4_RTS = log(f) - log(1-f);

load h10.mat
Z_10_AIS = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1), 100000)
% Z_10_SAMS = SAMS2(parameter_W, parameter_a, parameter_b, 1e-4)
% Z_10_TAP = TAP(parameter_W, parameter_b, parameter_a,1e-4)
% Z_10_RTS = RTS(parameter_W, parameter_a, parameter_bA_4_RTS, parameter_b, 0.25)

load h20.mat
Z_20_AIS = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1), 400000)
% Z_20_SAMS = SAMS(parameter_W, parameter_b, parameter_a)
% Z_20_TAP = TAP(parameter_W, parameter_b, parameter_a)
% Z_20_RTS = RTS(parameter_W, parameter_a, parameter_b, 0.3)

load h100.mat
Z_100_AIS = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1),500000)
% Z_100_SAMS = SAMS(parameter_W, parameter_b, parameter_a)
% Z_100_TAP = TAP(parameter_W, parameter_b, parameter_a)
% Z_100_RTS = RTS(parameter_W, parameter_a, parameter_b, 0.4)

load h500.mat
Z_500_AIS = AIS(parameter_W, parameter_a, parameter_a, parameter_b, parameter_b, ...
    zeros(size(parameter_W,1),1),1000000)
% Z_500_SAMS = SAMS(parameter_W, parameter_b, parameter_a)
% Z_500_TAP = TAP(parameter_W, parameter_b, parameter_a)
% Z_500_RTS = RTS(parameter_W, parameter_a, parameter_b, 0.5)


