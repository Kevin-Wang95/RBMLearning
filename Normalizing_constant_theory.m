clear all;close all;
load h10.mat

bin_add = @(x,y) dec2bin(bin2dec(x)+bin2dec(y));
z = 0;
v = '0';
h = '0';
v_aim = '1';
h_aim = '1';
for i = 1:length(parameter_b)
    v_aim = [v_aim '0'];
end
for j = 1:length(parameter_a)
    h_aim = [h_aim '0'];    
end
while(~strcmp(h ,h_aim))
    while(~strcmp(v ,v_aim))
        v_cal = v - '0';
        v_cal = [zeros(1,length(parameter_b)-length(v_cal)), v_cal];
        h_cal = h - '0';
        h_cal = [zeros(1,length(parameter_a)-length(h_cal)), h_cal];
        E = -v_cal*parameter_W*h_cal'-parameter_b*v_cal'-parameter_a*h_cal';
        z = z + exp(-E);
        v = bin_add(v, '1');
    end
    v = '0';
    h = bin_add(h, '1');
    disp(z);
end