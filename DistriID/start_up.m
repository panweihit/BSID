% startup script to make Matlab aware of the package
%
% (c) by Wei Pan, w.pan11@imperial.ac.uk, 2014 Mar 3
clear all;
clc;
close all;
%restoredefaultpath;
 
disp(['executing Wei Pan-BSID startup script...']);
me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located

addpath(mydir(1:end-1))

addpath([mydir,'demo'])
% addpath([mydir,'data'])
addpath([mydir,'solver'])
addpath([mydir,'subfunction'])
addpath([mydir,'inf'])
% if strcmp(computer, 'GLNXA64')
%     run('../cvx/cvx-a64/cvx_setup')
% elseif strcmp(computer, 'PCWIN64')
%     run('..\cvx\cvx-w64\cvx_setup')
% elseif strcmp(computer, 'MACI64')
%     run('../cvx/cvx-m64/cvx_setup')
% end

clear me mydir
