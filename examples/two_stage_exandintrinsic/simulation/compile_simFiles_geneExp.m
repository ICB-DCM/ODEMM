% This script generate the simulation files for the example of two stage
% gene expression including intrinsic and extrinsic noise. This requires
% the MATLAB Toolbox CERENA. The _syms.m files are generate using CERENA
% and the model definition files modelDef_geneExp.m for the reaction rate
% equations and the moment approximation, and modelDef_geneExp_extrinsic.m
% for the moment approximation with additional extrinsic noise.

clear all
close all
clc

% reaction rate equations
cvodewrap('geneExp_RRE','RRE_geneExp_RRE_syms')

% moment approximation
cvodewrap('geneExp_MA','MEC_2_LD_2_c_geneExp_MA_syms')

% moment approximation with additional extrinsic noise in all parameters
cvodewrap('geneExp_extrinsic','MEC_2_LD_2_c_geneExp_extrinsic_syms')


