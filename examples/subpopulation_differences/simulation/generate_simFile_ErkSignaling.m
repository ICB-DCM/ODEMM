% This script generates the functions for the NGF-induced Erk1/2 signaling
% using the toolbox AMICI.
clear all
close all
clc

amiwrap('ODEmodel_sPsET_loglog','ODEmodel_sPsET_loglog_syms',pwd)

