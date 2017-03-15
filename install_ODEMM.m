% Script that adds the required paths to the MATLAB search path.

ODEMM_path = fileparts(which('install_ODEMM.m'));
addpath(ODEMM_path);
addpath(genpath([ODEMM_path '/distributions/']));