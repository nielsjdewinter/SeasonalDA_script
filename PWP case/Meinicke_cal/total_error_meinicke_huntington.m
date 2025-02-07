function [err_T] = total_error_meinicke_huntington (err_D47,D47)
%This function propagates the D47 analytical error and the error of the D47-temperature calibration of
%Meinicke et al. (2020) [recalculated to the I-CDES scale in Meinicke et al. (2021)] on the final temperature

%The function needs two matlab files named slope_Niklas.mat and inter_Niklas.mat 
%These files contain the slopes and intercepts of the bootstrap York fit
%that was performed by Niklas Meinicke on the recalculated I-CDES calibration dataset of Meinicke et al. (2020)
%These files have to be in the same folder as the function 

%The math behind the approach was taken from the supplemental information of
%Huntington et al. (2009) RCMS

%%
%Need slope intercept pairs from bootstratp York Fit to calculate
%variacance-covariance matrix

load('slope_Niklas.mat');
load('inter_Niklas.mat');

slope=mean(slope_boot);
intercept=mean(inter_boot); 
CoefficientCovariance=cov(inter_boot,slope_boot);
var_intercept=CoefficientCovariance(1,1);
var_slope=CoefficientCovariance(2,2);
covari=CoefficientCovariance(1,2);

%% Calculate the intercept, slope and their variance and covariance
%convert analytical error in D47 to variance
err_D47=err_D47^2;

%% Calculate differentials required for variance-covariance matrix
d_slope = @(slope,D47,intercept) (500*sqrt(slope/(D47 - intercept)))/slope; %slope
 
d_47 = @(slope,D47,intercept)  -(500*(slope/(D47 - intercept))^(3/2))/slope; %D47
 
d_inter = @(slope,D47,intercept) (500*(slope/(D47 - intercept))^(3/2))/slope; %intercept
 
%% Establish variance covariance matrix
g = [d_inter(slope,D47,intercept) d_slope(slope,D47,intercept) d_47(slope,D47,intercept)];
 
v = [var_intercept,covari,0;
    covari,var_slope,0;
    0,0,err_D47;];

sigma=g*v*g'; 

%% The output of the function is the combined (total) error on temperature as 95% CI
err_T=sqrt(sigma)*2;
 
end

%References
%Meinicke, N. et al. A robust calibration of the clumped isotopes to temperature relationship for foraminifers. Geochim. Cosmochim. Acta 270, 160–183 (2020).
%Meinicke, N., Reimi, M. A., Ravelo, A. C. & Meckler, A. N. Coupled Mg/Ca and clumped isotope measurements indicate lack of substantial mixed layer cooling in the Western Pacific Warm Pool during the last∼ 5 million years. Paleoceanography and Paleoclimatology 36, e2020PA004115 (2021).
%Huntington, K. W. et al. Methods and limitations of ‘clumped’ CO2 isotope (Δ47) analysis by gas‐source isotope ratio mass spectrometry. Journal of Mass Spectrometry 44, 1318-1329 (2009).

