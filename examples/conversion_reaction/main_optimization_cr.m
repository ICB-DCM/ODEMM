% This is the main script for the first example, the conversion reaction.
% All Reaction Rate Equation (RRE) models are generated and their parameters
% are estimated by calling models_RRE.m. 
% All SP models are estimated, for which e.g. model_SP_k1.m indicates that k1 has
% an additional cell-to-cell variability modeled by the Sigma point
% approximation

clear all
close all
clc

models_RRE

model_SP_k3
model_SP_k2
model_SP_k1
model_SP_k1k3
model_SP_k1k2
model_SP_k2k3
model_SP_all

