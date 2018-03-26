% This script generates the functions for the conversion reaction using the
% toolbox AMICI

clear all
close all
clc

% generate simulation files used for sigma point approximation
amiwrap('CR_log','CR_log_syms',pwd)

% generate simulation files used for Reaction Rate Equations
amiwrap('CR','CR_syms',pwd)
