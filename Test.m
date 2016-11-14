close all
load test.mat
load h100.mat

z=zeros(8,1);

testbatchdata = permute(testbatchdata, [1,3,2]);
testbatchdata = reshape(testbatchdata,size(testbatchdata,1) ...
    *size(testbatchdata,2),size(testbatchdata,3));

f = (sum(testbatchdata,1)+5)/size(testbatchdata,1);
parameter_bA_4_RTS = log(f) - log(1-f);

for i = 1:20
   Z10(i) = RTS(parameter_W, parameter_a, parameter_bA_4_RTS, parameter_b,1.5)
end
