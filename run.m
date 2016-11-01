clear all;
load h10.mat
Z_theory = 0;

load test.mat
testbatchdata = reshape(permute(testbatchdata, [1,3,2]),size(testbatchdata,1) ...
    *size(testbatchdata,2),size(testbatchdata,3));

