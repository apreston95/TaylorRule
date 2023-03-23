clear all
clc

dynare NK
dynare NK_Uniform % This is the NK model with a Taylor rule and uniform priors on phi_pi to ensure no prior truncation when comparing models
dynare NK_RR