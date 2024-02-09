function [U, wopt] =  optDNA_nucleosome(S)

% Computes the minimum DNA nucleosome wrapping energy U and the 
% corresponding nucleosomal DNA configuration wopt
%
% Input: a DNA sequence of length 147, for example 
% S = 'ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT';

global whats
global stiff
global fra_with_const
global indc
global indw
global nc
global nn

indc = [ 5 15 25 36 46 56 66 76 87  97 107 117 127 138];
indw = [10 21 31 41 51 61 72 82 92 102 112 123 133 143];

nc = length([indc indw]);
nn = 19+72*25+nc/2;

addpath('cgDNA+')
params = load('Di-hemi-meth-hmeth.mat'); % cgDNA+ parameter set

[whats, stiff] = constructSeqParms(S, params);  

ref = load('initial_frames3.mat'); % Reference and initial configuration
fra_with_const = ref.bpl;
w_init = ref.w;

zvec  = z_from_w_bp_28ph_o_center(w_init); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Optimisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimoptions(@fminunc,'Algorithm','trust-region', 'SpecifyObjectiveGradient',true,'HessianFcn','objective','Display','off','MaxIter',1000,'TolFun',1e-16,'TolX',1e-16);
zeq  =  fminunc('nucleosome_dna_energy_grad_hess_bp_28ph_o_center',zvec,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Final: shape plus cgNA+energy %%%%%%%%%%%%%

[wopt,~]  = w_from_z_bp_28ph_o_center(zeq);
U = 0.5*(wopt - whats)'*stiff*(wopt - whats);

end











