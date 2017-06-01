% startup script to make Matlab aware of the package
% Copyright by Wei Pan, panweihit@gmail.com, 2017

clear all;
clc;
close all;


disp(['executing Wei Pan-Multiple Toobox startup script...']);
fprintf('You are using Bayesian System Identification toolbox for Multiple Dataset')
me = mfilename;                                            
mydir = which(me); mydir = mydir(1:end-2-numel(me));      

addpath(mydir(1:end-1)) 
addpath([mydir,'subfunction'])
addpath([mydir,'inf'])  
addpath([mydir,'demo'])
addpath([mydir,'solver'])

clear me mydir
