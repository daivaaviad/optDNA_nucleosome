function [eta,w,etapW,wpW,u,v,etapC,wpC] = vector2shapes(y)

%-------------------------------------------------------
% cgDNA function: [eta,w,etapW,wpW,u,v,etapC,wpC] = vector2shapes(y)
%-------------------------------------------------------
% This function re-orders the ground-state coordinates.
%
% Input:
%
%    y    overall coordinate vector
%         [size N x 1].
%
%
% Output:
%
%    eta  list of intra-basepair rotational coords
%         (Buckle,Propeller,Opening) along molecule
%         [size nbp x 3]
%
%    w    list of intra-basepair translational coords
%         (Shear,Stretch,Stagger) along molecule
%         [size nbp x 3]
%
%    u    list of inter-basepair rotational coords
%         (Tilt,Roll,Twist) along molecule
%         [size (nbp-1) x 3]
%
%    v    list of inter-basepair translational coords
%         (Shift,Slide,Rise) along molecule
%         
%    etapW, etapC  list of Watson strand and Crick strand phosphate
%         rotational coords along molecule
%         [size (nbp-1) x 3] 
%       
%    wpW, wpC  list of Watson strand and Crick strand phosphate
%         translational coords along molecule
%         [size (nbp-1) x 3] 
%
%
%    where N = 24*nbp - 18 and nbp is the length
%    of the DNA sequence (number of basepairs).
%
%
% If you find this code useful, please cite:
%
% Patelli, A. S., 2019. A sequence-dependent coarse-grain model 
% of B-DNA with explicit description of bases and phosphate
% groups parametrised from large scale Molecular Dynamics simulations. 
% Ph.D. thesis, EPFL.
%
%
% D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks.
%  cgDNA: a software package for the prediction of sequence-dependent
%  coarse-grain free energies of B-form DNA. Submitted (2014).
%
%-------------------------------------------------------

N   = numel(y);
nbp = (N+18)/24;

q	= reshape(y', 6, 4*nbp-3)';
intra = q(1:4:end, :);  % nbp
pho_C = q(2:4:end, :);  % nbp - 1
inter = q(3:4:end, :);  % nbp - 1
pho_W = q(4:4:end, :);  % nbp - 1

eta = intra(:, 1:3);
w  = intra(:, 4:6);

etapW = pho_W(:,1:3);
wpW   = pho_W(:,4:6);

u  = inter(:, 1:3);
v = inter(:, 4:6);

etapC = pho_C(:,1:3);
wpC   = pho_C(:,4:6);

end
