%This script will calculate the total error on the three D47-based temperatures 
%(background, slope, and hyperthermal peak) including both analytical and calibration uncertainties 

%Type name of comma delimited excel file here
%Inputs are the mean D47 values of the three bins (first column) and
%the D47 analytical error of the three bins as 68% CI (second column)
input_file='test.csv';
data = csvread(input_file);
 
D47 = (data(:,1));
err_D47 = (data(:,2));
 
for i=1:length(data)
    %Function to propagate analytical and calibration uncertainties on the
    %three D47-based temperatures
    %Function returns combined (total) error on the three temperatures as 95% CI 
    err_T(i) = total_error_meinicke_huntington (err_D47(i),D47(i))
    
end

 