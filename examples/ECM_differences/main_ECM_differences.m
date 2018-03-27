% Main script for the analysis of the influence of extracellular scaffolds
% on pain signaling in sensory neurons.

clear all
close all
clc

% Generate the hierarchical population models for all 128 models
generate_ECM_models

% Run the fittings for all 128 models
run_fittings_ECM

% Run the profile likelhiood analysis for the model accounting for
% differenes in TrkA activity, Erk levels and Erk inactivation
run_profiles_ECM

