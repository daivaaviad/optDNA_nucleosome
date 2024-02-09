function [shapes, stiff] = constructSeqParms(seq, paramset)
% -------------------------------------------------------------------------
% cgDNAp function: [shapes,stiff] = constructSeqParms(seq, paramset)
% -------------------------------------------------------------------------
% This function constructs the ground-state coordinate
% vector and stiffness matrix in non-dimensional Curves+
% form for a given sequence, using the specified parameter
% set in params.  
%
%
% Input: 
%
%   seq     sequence along reference strand;
%
%   params  a cgDNAplus parameter set. See Note 1 for details. 
%
%
% Output:
%
%   shapes  ground-state coordinate vector 
%           [size N x 1]
%
%   stiff   ground-state stiffness matrix
%           [size N x N]
%
%   where N = 24*nbp - 18 and nbp is the length 
%   of the sequence seq (number of basepairs). 
%
%
% Note 1: The cgDNAplus paramset is given as a MATLAB structure which
%         contains the following fields : 
%
%         sigma_int  : Contains the values of sigma for the iterior dimers
%         sigma_end5 : Contains the values of sigma for the 5' end dimers
%         sigma_end3 : Contains the values of sigma for the 3' end dimers
%
%         stiff_int  : Contains the stiff. blocks for the iterior dimers
%         stiff_end5 : Contains the stiff. blocks for the 5' end dimers
%         stiff_end3 : Contains the stiff. blocks for the 3' end dimers
%         
%         Each field contains itself one field corresponding to each dimer.
%
%
% If you find this code useful, please cite:
%
%  Patelli, A. S., 2019. A sequence-dependent coarse-grain model 
%  of B-DNA with explicit description of bases and phosphate
%  groups parametrised from large scale Molecular Dynamics simulations.
%  Ph.D. thesis, EPFL.
%
%  Petkeviciute, D., M. Pasi, O. Gonzalez, and J. Maddocks, 2014. 
%  cgDNA: a software package for the prediction of sequence-dependent 
%  coarse-grain free energies of B-form DNA. Nucleic Acids Res. 42:e153–e153.
%
% -------------------------------------------------------------------------

% Initialize variables
seq = upper(seq);
nbp = numel(seq);
N =24*nbp-18 ;
stiff_cell = cell(nbp-1,1);
sigma = zeros(N,1);

dimer = repmat(seq, [nbp-1 1]) ; 
dimer = [diag(dimer),diag(dimer,1)] ;

%5' endblocks
[I,J,stiff_block] = find(paramset.stiff_end5.(dimer(1,:)) );
stiff_cell{1}= [I,J,stiff_block];
sigma(1:36) = paramset.sigma_end5.(dimer(1,:))  ;


%Stiffness Matrix and sigma vector's core
for i = 2:nbp-2
    k = (i-2)*24+18;
    
    %Sparsity pattern indices and corresponding dimer block values are
    %stored
    [I,J,stiff_block] =  find( paramset.stiff_int.(dimer(i,:)));
    stiff_cell{i}= [I+ k*ones(size(I)),J+ k*ones(size(J)),stiff_block];
    sigma(k+1:k+42,1) = sigma(k+1:k+42,1) + paramset.sigma_int.(dimer(i,:));
    
end

k = (nbp-3)*24 + 18;

%3' endblocks
[I,J,stiff_block] = find(paramset.stiff_end3.(dimer(nbp-1,:)));
stiff_cell{nbp-1}= [I + k*ones(size(I)),J + k*ones(size(J)),stiff_block];

sigma(k+1:k+36,1) = sigma(k+1:k+36,1)+ paramset.sigma_end3.(dimer(nbp-1,:)) ;


%Assembly of the Sparse Stiffness matrix
IJV= cell2mat(stiff_cell);
stiff = sparse(IJV(:,1),IJV(:,2),IJV(:,3),N,N);

shapes = stiff\sigma;

end