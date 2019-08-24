%Analysis of brownian motion simulation data in generated files
tic
Number_of_files = 60;
filename_pre = 'C:\Users\d-w-h\Desktop\Home\Brownian Motion Simulation results\BrownianMotionSim';
filename_pos = 'Data.xlsx';

%Read one file for data retrieval
index = num2str(1);
file_name_front = strcat(filename_pre, index);
file_name = strcat(file_name_front, filename_pos);
P = xlsread(file_name,2,'A1:C20002');

Length_data = length(P(:,1));

%Read all files
dR_sq_avg_overall = zeros(Length_data,1);
for i = 1:Number_of_files
    index = num2str(i);
    file_name_front = strcat(filename_pre, index);
    file_name = strcat(file_name_front, filename_pos)
    M = xlsread(file_name,2,'A1:C20002');

    t = M(:,1);
    dR_sq_avg = M(:,2);
    D_avg = M(:,3);
    
    dR_sq_avg_overall = dR_sq_avg + dR_sq_avg_overall;
end

%Calculate mean squared distance
dR_sq_avg_overall = dR_sq_avg_overall / Number_of_files;

%Calculating average diffusion coefficient
T = P(:,1);
D_avg_overall = dR_sq_avg_overall(end) / 4 / T(end);

%Expected mean squared distance based on diffusion coefficient
dR_sq_avg_exp = 4 * D_avg_overall * T;

%Plot figures
hold on
plot(T, dR_sq_avg_overall)
plot(T, dR_sq_avg_exp)
hold off

%Export data
A = [T, dR_sq_avg_overall];
filename = 'C:\Users\d-w-h\Desktop\Home\Brownian Motion Simulation results\BrownianMotionAVGData.xlsx';
xlswrite(filename,A,1,'A1:B25000')

toc