function y = shapes2vector(eta,w,etapW,wpW,u,v,etapC,wpC)

%-------------------------------------------------------
% cgDNA function: y = shapes2vector(eta,w,etapW,wpW,u,v,etapC,wpC)
%-------------------------------------------------------
% This function re-orders the ground-state coordinates.
%
% Input: 
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
%         [size (nbp-1) x 3].
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
% Output:
%
%    y    overall coordinate vector 
%         [size N x 1]
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
% cgDNA: a software package for the prediction of sequence-dependent
% coarse-grain free energies of B-form DNA. Submitted (2014).
%
%-------------------------------------------------------

    nbp = size(eta,1);
    N = 24*nbp - 18;
    
    rotat = [ eta ; etapC ; u ; etapW ];
    trasl = [   w ; wpC ; v ; wpW ];
    q = [ rotat, trasl ];
    indices = repmat( 1:nbp-1 , [4 1] ) + repmat( [ 0 ; nbp ; 2*(nbp-1)+1 ; 3*(nbp-1)+1] , [1 nbp-1] ) ;
    indices = [indices(:); nbp];
    
    y = reshape( q(indices,:)', 1, N )';
    
end
