function [energy,grad,hess] = nucleosome_dna_energy_grad_hess_bp_28ph_o_center(zvec)

% Computes the cgNA+ energy with an added penalty function on 
% phosphate positions for wrapping DNA into a nucleosome, as well
% as the corresponding gradient and Hessian.
%
% Input: zvec - a vector of hybrid coordinates


%% ----------- Initialisation -------------------

global whats
global stiff
global fra_with_const
global indc
global indw
global nc
global nn

nuc = 1;
ph_penalty_weight = 0.01;
qnorm_penalty_weight = 1000;
wi =  [2 2 2 3 5 7 10 10 7 5 3 2 2 2];

zlen = length(zvec)+7;
wlen = length(whats);

[wvec, zel] = w_from_z_bp_28ph_o_center(zvec);

qbp = zel.qbp;
obp = zel.obp;
qbc = zel.qbc;
obc = zel.obc;
qbw = zel.qbw;
obw = zel.obw;
opw = zel.opw;
opc = zel.opc;

nbp = size(qbp,1);

% Constant matrices
b  = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3)  = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1)  = -1; b(3,3,4) = 1; b(3,4,3) = -1;

B1 = reshape(b(1,:,:),4,4);

[i33, j33, i44, j44, i34, j34, i43, j43] = index_matrices();

if nuc
  const = fra_with_const;
end    

%% ----------- Energy + penalty -----------------

energy = (wvec - whats)'*stiff*(wvec - whats)/2;

% Nucleosome penalty

if nuc %set nuc = 0 for testing
%indc =  [6 16 25 36 46 56 66 76 87 97 108 118 127 137];
%indw =  [11 21 30 40 51 61 72 82 92 102 112 123 132 142];

for i = indc % constraints on positions of C phosphates
    pos_variable   = opc(i,:)'; 
    pos_constraint = const(i).rpc;
    energy = energy + ph_penalty_weight*wi(i == indc)*(norm(pos_variable - pos_constraint)^2);
end

for i = indw % constraints on positions of W phosphates
    pos_variable   = opw(i,:)'; 
    pos_constraint = const(i).rpw;
    energy = energy + ph_penalty_weight*wi(i == indw)*(norm(pos_variable - pos_constraint)^2);  
end
end

% Quaternion norm penalty    
for i=1:nbp 
    energy = energy + qnorm_penalty_weight * ((qbp(i,:)*qbp(i,:)'-1)^2);  
    if sum(i == indc) 
         energy = energy + qnorm_penalty_weight * ((qbc(i,:)*qbc(i,:)'-1)^2);      
    elseif sum(i == indw) 
         energy = energy + qnorm_penalty_weight * ((qbw(i,:)*qbw(i,:)'-1)^2);     
    end    
end
                                  
%% ----------- Gradient -------------------------

% First compute the Jacobian dw/dz
nc = length([indc indw]);
N  = (nbp-1)*22+(nbp-nc-1)*18+6+nc*37; 
I  = zeros(3,N);  % Jacobian row indices
K  = zeros(3,N);  % Jacobian column indices
Jb = zeros(3,N);  % Jacobian blocks
j = 1;

jb = eye(3);

% Inter derivatives
k2 = 0;
for i = 1:nbp-1
    
     k1 = (i-1)*24; 
     if i==1
         step = 19;
     elseif sum(i==[indc indw])
         step = 26;
     else
         step = 25;
     end
     
     [jb1, jb2] = dcay_dq(qbp(i,:)', qbp(i+1,:)', wvec(k1+13:k1+15), b);
     I(:,j:j+7)  = [k1+13+i34, k1+13+i34];
     K(:,j:j+7)  = [k2+4+j34,  k2+4+step+j34];
     Jb(:,j:j+7) = [jb1, jb2]; 
     j = j+8;
       
     [jb3, jb4] = dtrans_doq(qbp(i,:)', qbp(i+1,:)', obp(i,:)', obp(i+1,:)',b, 0);   
     I(:,j:j+13)  = [k1+16+i33, k1+16+i34, k1+16+i33, k1+16+i34];
     K(:,j:j+13)  = [k2+1+j33,  k2+4+j34,  k2+1+step+j33,  k2+4+step+j34];
     Jb(:,j:j+13) = [-jb3,  jb4,  jb3,  jb4]; 
     j = j+14;   
     
     k2 = k2+step;
end    

k2 = 0;
for i = 1:nbp
    
  k1 = (i-1)*24;  
  
  if sum(i==indc) % Constrained C phosphates   
      
      % Intra derivatives
      
      cay = cayley_from_q(qbc(i,:)', qbp(i,:)', b);    
      [jb1, jb2] = dcay_dq(qbc(i,:)', qbp(i,:)', cay, b);
      
      dc = 200/(100 - cay'*cay) * eye(3) + 400/(100 - cay'*cay)^2 * (cay*cay');
      djb1 = dc * jb1;
      djb2 = dc * jb2;
      
      I(:,j:j+7)  = [k1+1+i34, k1+1+i34];
      K(:,j:j+7)  = [k2+11+j34, k2+4+j34];
      Jb(:,j:j+7) = [djb1, djb2];
      j = j+8;
      
      [jb3, jb4] = dtrans_doq(qbp(i,:)', qbp(i,:)', obc(i,:)', obp(i,:)', b, 1);
      I(:,j:j+9)  = [k1+4+i33, k1+4+i34, k1+4+i33];
      K(:,j:j+9)  = [k2+1+j33, k2+4+j34, k2+8+j33];
      Jb(:,j:j+9) = [2*jb3, 2*jb4, -2*jb3];
      j=j+10;
      
      % C phos derivatives wrt obc, qbc, opc 
      % Rotations (identity matrices)
      I(:,j:j+2)  = k1+7+i33;
      K(:,j:j+2)  = k2+18+j33;
      Jb(:,j:j+2) = jb; 
      j = j+3;
  
      % Translations (B1*qbc is flipping the base frame)
      [jb3, jb4] = dtrans_doq(B1*qbc(i,:)', qbc(i,:)', obc(i,:)', opc(i,:)', b, 1); %second argument unused for phosphates
      I(:,j:j+9)  = [k1+10+i33, k1+10+i34, k1+10+i33];
      K(:,j:j+9)  = [k2+8+j33,  k2+11+j34, k2+15+j33];
      Jb(:,j:j+9) = [-jb3, jb4*B1, jb3]; 
      j=j+10;
      
      % W phos derivatives (identity matrices)
      I(:,j:j+5)  = [k1-5+i33, k1-2+i33];
      K(:,j:j+5)  = [k2+21+j33, k2+24+j33];
      Jb(:,j:j+5) = [jb, jb];
      j = j+6;
      
      k2 = k2+26;
      
  elseif sum(i==indw) % Constrained W phosphates 
      
      % Intra derivatives
      cay = cayley_from_q(qbp(i,:)', qbw(i,:)', b);
      [jb1, jb2] = dcay_dq(qbp(i,:)', qbw(i,:)', cay, b);
      
      dc = 200/(100 - cay'*cay) * eye(3) + 400/(100 - cay'*cay)^2 * (cay*cay');
      djb1 = dc * jb1;
      djb2 = dc * jb2;
      
      I(:,j:j+7)  = [k1+1+i34, k1+1+i34];
      K(:,j:j+7)  = [k2+4+j34, k2+11+j34];
      Jb(:,j:j+7) = [djb1, djb2];
      j = j+8;
  
      [jb3, jb4] = dtrans_doq(qbp(i,:)', qbp(i,:)', obp(i,:)', obw(i,:)', b, 1);
      I(:,j:j+9)  = [k1+4+i33, k1+4+i34, k1+4+i33];
      K(:,j:j+9)  = [k2+1+j33, k2+4+j34, k2+8+j33];
      Jb(:,j:j+9) = [-2*jb3, 2*jb4, 2*jb3];
      j=j+10;
      
      % W phos derivatives wrt obw, qbw, opw
      % Rotations (identity matrices)
      I(:,j:j+2)  = k1-5+i33;
      K(:,j:j+2)  = k2+18+j33;
      Jb(:,j:j+2) = jb; 
      j = j+3;
  
      % Translations
      [jb3, jb4] = dtrans_doq(qbw(i,:)', qbw(i,:)', obw(i,:)', opw(i,:)',b, 1);
      I(:,j:j+9)  = [k1-2+i33, k1-2+i34, k1-2+i33];
      K(:,j:j+9)  = [k2+8+j33,  k2+11+j34, k2+15+j33];
      Jb(:,j:j+9) = [-jb3, jb4, jb3]; 
      j=j+10;
      
      % C phos derivatives (identity matrices)
      I(:,j:j+5)  = [k1+7+i33, k1+10+i33];
      K(:,j:j+5)  = [k2+21+j33, k2+24+j33];
      Jb(:,j:j+5) = [jb, jb];
      j = j+6;
      k2 = k2+26;
      
  else % Unconstrained part   
  
     % Intra derivatives (identity matrices)
     I(:,j:j+5)  = [k1+1+i33, k1+4+i33];
     K(:,j:j+5)  = [k2+8+j33, k2+11+j33];
     Jb(:,j:j+5) = [jb, jb];
     j = j+6;
    
     if i<nbp
        % C phos derivatives (identity matrices)
        I(:,j:j+5)  = [k1+7+i33, k1+10+i33];
        K(:,j:j+5)  = [k2+14+j33, k2+17+j33];
        Jb(:,j:j+5) = [jb, jb];
        j = j+6;  
     end
     
     if i>1
         if i<nbp
           % W phos derivatives (identity matrices)
           I(:,j:j+5)  = [k1+19-24+i33, k1+22-24+i33];
           K(:,j:j+5)  = [k2+20+j33, k2+23+j33];
           Jb(:,j:j+5) = [jb, jb];
           j = j+6;
         else
           I(:,j:j+5)  = [k1-5+i33, k1-2+i33];
           K(:,j:j+5)  = [k2+14+j33, k2+17+j33];
           Jb(:,j:j+5) = [jb, jb];
         end    
     end  
   
     if i==1
        k2 = 19;
     else   
        k2 = k2+25;
     end   
   
  end
end

% Compute the gradient
J = sparse(I, K, Jb, wlen, zlen);
grad0 = stiff*(wvec - whats);
grad = J'*grad0;

% Add constraints

if nuc
% Nucleosome penalty
k = 19;
for i = 2:nbp
    if sum(i == indc) 
      pos_variable   = opc(i,:)';
      pos_constraint = const(i).rpc;
      grad(k+15:k+17) = grad(k+15:k+17) ...
                       + 2*ph_penalty_weight*wi(i == indc)*(pos_variable - pos_constraint); 
      k=k+26;
      
    elseif sum(i == indw) 
      pos_variable   = opw(i,:)';
      pos_constraint = const(i).rpw;
      grad(k+15:k+17) = grad(k+15:k+17) ...
                        + 2*ph_penalty_weight*wi(i == indw)*(pos_variable - pos_constraint); 
      k=k+26;
      
    else
      k=k+25;
    end  
end
end

% Quaternion norm penalty

grad(4:7) = grad(4:7) + 4*qnorm_penalty_weight * (qbp(1,:)*qbp(1,:)'-1)*qbp(1,:)'; 
k = 19;
for i=2:nbp 
    grad(k+4:k+7) = grad(k+4:k+7) + 4*qnorm_penalty_weight * (qbp(i,:)*qbp(i,:)'-1)*qbp(i,:)';  
    if sum(i == indc) 
        grad(k+11:k+14) = grad(k+11:k+14) + 4*qnorm_penalty_weight * (qbc(i,:)*qbc(i,:)'-1)*qbc(i,:)';
        k = k+26;    
    elseif sum(i == indw) 
        grad(k+11:k+14) = grad(k+11:k+14) + 4*qnorm_penalty_weight * (qbw(i,:)*qbw(i,:)'-1)*qbw(i,:)';
        k = k+26; 
    else
        k = k+25; 
    end    
end

grad = grad([1:nn,nn+8:end]); %the base pair number 74 (the dyad point) is fixed
 
%% ----------- Hessian --------------------------

% First compute the derivatives d^2 w_k1 / d(zvec_k2)^2 
%   (a list of sparse Hessian matrices with very few nonzero elements)
hw_list  = struct('I', [],'K',[], 'H', cell(wlen,1));

ind = zeros(10,1);
iind = 0;

k2 = 0;
for i=1:nbp-1
    
     k1 = (i-1)*24;   
     if i==1
        step = 19;
     elseif sum(i==[indc indw])
        step = 26;
     else
        step = 25;
     end
                  
     % Inter derivatives:

     [hb11, hb22, hb12] = ddcay_dqdq(qbp(i,:)', qbp(i+1,:)', wvec(k1+13:k1+15), b);
     for j = 1:3  
          
        hb11j = reshape(hb11(j,:,:),4,4);
        hb22j = reshape(hb22(j,:,:),4,4);
        hb12j = reshape(hb12(j,:,:),4,4)';
            
        hw_list(k1+12+j).I = [k2+4+i44, k2+4+step+i44, k2+4+i44,  k2+4+step+j44];                    
        hw_list(k1+12+j).K = [k2+4+j44, k2+4+step+j44, k2+4+step+j44, k2+4+i44];                    
        hw_list(k1+12+j).H = [hb11j, hb22j, hb12j, hb12j];  
                             
     end 
                
     [hbo, hbqq] = ddtrans_ddoq(qbp(i,:)', qbp(i+1,:)', obp(i,:)', obp(i+1,:)', b, 0);
     for j = 1:3 
          
          hboj = reshape(hbo(j,:,:),4,3);
          hbqqj = reshape(hbqq(j,:,:),4,4);
          
          hw_list(k1+15+j).I = [k2+4+i44,  k2+4+step+i44, k2+4+i44,  k2+4+step+j44, ...                     
                                    k2+4+i43,  k2+4+step+i43, k2+4+i43, k2+4+step+i43, ...
                                    k2+1+j43,  k2+1+j43,  k2+1+step+j43,  k2+1+step+j43];
          hw_list(k1+15+j).K = [k2+4+j44,  k2+4+step+j44, k2+4+step+j44, k2+4+i44, ...
                                    k2+1+j43,  k2+1+j43,  k2+1+step+j43,  k2+1+step+j43, ...
                                    k2+4+i43,  k2+4+step+i43, k2+4+i43, k2+4+step+i43];
          hw_list(k1+15+j).H = [hbqqj, hbqqj, hbqqj, hbqqj, ...  
                                    -hboj, -hboj, hboj, hboj, ...
                                    -hboj, -hboj, hboj, hboj];      
     end  
     
     ind(iind+1:iind+6) = k1+13:k1+18;
     iind = iind+6;
     
     % Only constrained intra and phosphate second derivatives are nonzero
     
     if sum(i==indc) % Constrained C phosphates   
      
       % Intra derivatives
       
       cay = cayley_from_q(qbc(i,:)', qbp(i,:)', b); 
       [hb11, hb22, hb12] = ddcay_dqdq(qbc(i,:)', qbp(i,:)', cay, b);
       
       [jb1, jb2] = dcay_dq(qbc(i,:)', qbp(i,:)', cay, b);
       cc = 100 - cay'*cay;
       dsum11 = (cay(1)*reshape(hb11(1,:,:),4,4) + ...
                      cay(2)*reshape(hb11(2,:,:),4,4) + cay(3)*reshape(hb11(3,:,:),4,4));
       dsum22 = (cay(1)*reshape(hb22(1,:,:),4,4) + ...
                      cay(2)*reshape(hb22(2,:,:),4,4) + cay(3)*reshape(hb22(3,:,:),4,4));  
       dsum12 = (cay(1)*reshape(hb12(1,:,:),4,4) + ...
                      cay(2)*reshape(hb12(2,:,:),4,4) + cay(3)*reshape(hb12(3,:,:),4,4));    
       
       for j = 1:3 
           
         hb11j = reshape(hb11(j,:,:),4,4);
         hb22j = reshape(hb22(j,:,:),4,4);
         hb12j = reshape(hb12(j,:,:),4,4)';           
                  
         dd = zeros(3);
         dd(:,j) = cay;
         dd = dd + cay(j)*eye(3);
           
         hb11jc =  400*cc^(-2) * jb1(j,:)'*cay'*jb1 + 200/cc * hb11j ...
                + 1600*cc^(-3) * cay(j)*jb1'*(cay*cay')*jb1 ...
                +  400*cc^(-2) * (jb1'*dd*jb1 + cay(j)*dsum11);
            
         hb22jc =  400*cc^(-2) * jb2(j,:)'*cay'*jb2 + 200/cc * hb22j ...
                + 1600*cc^(-3) * cay(j)*jb2'*(cay*cay')*jb2 ...
                +  400*cc^(-2) * (jb2'*dd*jb2 + cay(j)*dsum22);   
          
         hb12jc =  400*cc^(-2) * jb1(j,:)'*cay'*jb2 + 200/cc * hb12j ...
                + 1600*cc^(-3) * cay(j)*jb1'*(cay*cay')*jb2 ...
                +  400*cc^(-2) * (jb1'*dd*jb2 + cay(j)*dsum12'); 

         hw_list(k1+j).I = [k2+11+i44, k2+4+i44, k2+11+i44,  k2+4+j44];
         hw_list(k1+j).K = [k2+11+j44, k2+4+j44, k2+4+j44, k2+11+i44];
         hw_list(k1+j).H = [hb11jc, hb22jc, hb12jc, hb12jc];   
       end    
  
       [hbo, hbqq] = ddtrans_ddoq(qbp(i,:)', qbp(i,:)', obc(i,:)', obp(i,:)', b, 1);
       for j = 1:3  
          %(d^2 trans / dodo = 0;)
          hboj = reshape(hbo(j,:,:),4,3);
          hbqqj = reshape(hbqq(j,:,:),4,4);
          hw_list(k1+3+j).I = [k2+4+i44, k2+4+i43,  k2+4+i43,  ...
                                k2+1+j43,  k2+8+j43];
          hw_list(k1+3+j).K = [k2+4+j44, k2+1+j43,  k2+8+j43,  ...
                                k2+4+i43,  k2+4+i43];
          hw_list(k1+3+j).H = [2*hbqqj, 2*hboj,  -2*hboj,  ...
                                2*hboj, -2*hboj];                       
       end 
      
       % C phos derivatives wrt obc, qbc, opc
       % Rotations are zero
        
       % Translations
       [hbo, hbqq] = ddtrans_ddoq(B1*qbc(i,:)', qbc(i,:)', obc(i,:)', opc(i,:)',b, 1);
       for j = 1:3  
           hboj = reshape(hbo(j,:,:),4,3);
           hbqqj = reshape(hbqq(j,:,:),4,4);
           hw_list(k1+9+j).I = [k2+11+i44, k2+11+i43,  k2+11+i43,  ...
                                          k2+8+j43,  k2+15+j43];
           hw_list(k1+9+j).K = [k2+11+j44, k2+8+j43,  k2+15+j43,  ...
                                          k2+11+i43,  k2+11+i43, ];
           hw_list(k1+9+j).H = [B1'*hbqqj*B1, -B1'*hboj,  B1'*hboj, ...
                                          -B1'*hboj, B1'*hboj];                  
       end 
       ind(iind+1:iind+12) = k1+1:k1+12;
       iind = iind+12;
      
     elseif sum(i==indw) % Constrained W phosphates   
      
      % Intra derivatives
      
       cay = cayley_from_q(qbp(i,:)', qbw(i,:)', b); 
       [hb11, hb22, hb12] = ddcay_dqdq(qbp(i,:)', qbw(i,:)', cay, b);
       
       [jb1, jb2] = dcay_dq(qbp(i,:)', qbw(i,:)', cay, b);
       cc = 100 - cay'*cay;
       dsum11 = (cay(1)*reshape(hb11(1,:,:),4,4) + ...
                      cay(2)*reshape(hb11(2,:,:),4,4) + cay(3)*reshape(hb11(3,:,:),4,4));
       dsum22 = (cay(1)*reshape(hb22(1,:,:),4,4) + ...
                      cay(2)*reshape(hb22(2,:,:),4,4) + cay(3)*reshape(hb22(3,:,:),4,4));  
       dsum12 = (cay(1)*reshape(hb12(1,:,:),4,4) + ...
                      cay(2)*reshape(hb12(2,:,:),4,4) + cay(3)*reshape(hb12(3,:,:),4,4));    
       
       for j = 1:3 
           
         hb11j = reshape(hb11(j,:,:),4,4);
         hb22j = reshape(hb22(j,:,:),4,4);
         hb12j = reshape(hb12(j,:,:),4,4)';           
                  
         dd = zeros(3);
         dd(:,j) = cay;
         dd = dd + cay(j)*eye(3);
           
         hb11jc =  400*cc^(-2) * jb1(j,:)'*cay'*jb1 + 200/cc * hb11j ...
                + 1600*cc^(-3) * cay(j)*jb1'*(cay*cay')*jb1 ...
                +  400*cc^(-2) * (jb1'*dd*jb1 + cay(j)*dsum11);
            
         hb22jc =  400*cc^(-2) * jb2(j,:)'*cay'*jb2 + 200/cc * hb22j ...
                + 1600*cc^(-3) * cay(j)*jb2'*(cay*cay')*jb2 ...
                +  400*cc^(-2) * (jb2'*dd*jb2 + cay(j)*dsum22);   
          
         hb12jc =  400*cc^(-2) * jb1(j,:)'*cay'*jb2 + 200/cc * hb12j ...
                + 1600*cc^(-3) * cay(j)*jb1'*(cay*cay')*jb2 ...
                +  400*cc^(-2) * (jb1'*dd*jb2 + cay(j)*dsum12'); 

         hw_list(k1+j).I = [k2+4+i44, k2+11+i44, k2+4+i44,  k2+11+j44];
         hw_list(k1+j).K = [k2+4+j44, k2+11+j44, k2+11+j44, k2+4+i44];
         hw_list(k1+j).H = [hb11jc, hb22jc, hb12jc, hb12jc];   
       end    
        
       [hbo, hbqq] = ddtrans_ddoq(qbp(i,:)', qbp(i,:)', obp(i,:)', obw(i,:)', b, 1);
       for j = 1:3  
         %(d^2 trans / dodo = 0;)
         hboj = reshape(hbo(j,:,:),4,3);
         hbqqj = reshape(hbqq(j,:,:),4,4);
         hw_list(k1+3+j).I = [k2+4+i44, k2+4+i43,  k2+4+i43,  ...
                               k2+1+j43,  k2+8+j43];
         hw_list(k1+3+j).K = [k2+4+j44, k2+1+j43,  k2+8+j43,  ...
                               k2+4+i43,  k2+4+i43];
         hw_list(k1+3+j).H = [2*hbqqj, -2*hboj,  2*hboj,  ...
                               -2*hboj, 2*hboj];                       
       end 
      
       % W phos derivatives wrt obw, qbw, opw    
       % Rotations are zero
        
       % Translations
          
       [hbo, hbqq] = ddtrans_ddoq(qbw(i,:)', qbw(i,:)', obw(i,:)', opw(i,:)', b, 1);
       for j = 1:3  
          hboj = reshape(hbo(j,:,:),4,3);
          hbqqj = reshape(hbqq(j,:,:),4,4);
          hw_list(k1+21+j-24).I = [k2+11+i44, k2+11+i43, k2+11+i43,  ...
                                            k2+8+j43,  k2+15+j43];
          hw_list(k1+21+j-24).K = [k2+11+j44, k2+8+j43,  k2+15+j43, ...
                                            k2+11+i43,  k2+11+i43];
          hw_list(k1+21+j-24).H = [hbqqj, -hboj, hboj, ...
                                            -hboj, hboj];                 
       end 
       ind(iind+1:iind+6) = k1+1:k1+6;
       ind(iind+7:iind+12) = k1-5:k1;
       iind = iind+12;
         
     end  
     k2 = k2+step;
      
end     

% Hessian of the cgDNA+ energy wrt zvec

hess = spalloc(zlen, zlen, 196620);
for k = ind'
     hess = hess + grad0(k) * sparse(hw_list(k).I, hw_list(k).K, hw_list(k).H, zlen, zlen); 
end

hess =  J' * stiff * J + hess;


% Nucleosome penalty
if nuc
k = 19;
for i = 2:nbp
    if sum(i == indc) 
      hess(k+15:k+17,k+15:k+17) = hess(k+15:k+17,k+15:k+17) ...
                                + 2*ph_penalty_weight*wi(i == indc)*eye(3); 
      k=k+26;
      
    elseif sum(i == indw) 
      hess(k+15:k+17,k+15:k+17) = hess(k+15:k+17,k+15:k+17) ...
                                + 2*ph_penalty_weight*wi(i == indw)*eye(3); 
      k=k+26;
      
    else
      k=k+25;
    end  
end
end


% Quaternion norm penalty
hess(4:7,4:7) = hess(4:7,4:7) + 4*qnorm_penalty_weight * ...   
                                       ((qbp(1,:)*qbp(1,:)'-1)*eye(4) + 2*qbp(1,:)'*qbp(1,:)); 
k = 19;
for i=2:nbp 
    hess(k+4:k+7,k+4:k+7) = hess(k+4:k+7,k+4:k+7) + 4*qnorm_penalty_weight * ...
                                       ((qbp(i,:)*qbp(i,:)'-1)*eye(4) + 2*qbp(i,:)'*qbp(i,:)); 
    if sum(i == indc) 
        hess(k+11:k+14,k+11:k+14) = hess(k+11:k+14,k+11:k+14) + 4*qnorm_penalty_weight * ...
                                       ((qbc(i,:)*qbc(i,:)'-1)*eye(4) + 2*qbc(i,:)'*qbc(i,:)); 
        k = k+26;    
    elseif sum(i == indw) 
        hess(k+11:k+14,k+11:k+14) = hess(k+11:k+14,k+11:k+14) + 4*qnorm_penalty_weight * ...
                                       ((qbw(i,:)*qbw(i,:)'-1)*eye(4) + 2*qbw(i,:)'*qbw(i,:)); 
        k = k+26; 
    else
        k = k+25; 
    end    
end

hess = hess([1:nn,nn+8:end], [1:nn,nn+8:end]); %the base pair number 74 (the dyad point) is fixed


end

%% ----------- Various functions ----------------

function [i33, j33, i44, j44, i34, j34, i43, j43] = index_matrices()

% Indices i and j of entries inside 3x3, 4x4, 3x4 and 4x4 blocks 
% for producing sparse matrices

  i33 = [0 0 0; 1 1 1; 2 2 2];
  j33 = [0 1 2; 0 1 2; 0 1 2]; 
  
  i44 = [0 0 0 0; 1 1 1 1; 2 2 2 2; 3 3 3 3];
  j44 = [0 1 2 3; 0 1 2 3; 0 1 2 3; 0 1 2 3]; 
  
  i34 = [0 0 0 0; 1 1 1 1; 2 2 2 2];
  j34 = [0 1 2 3; 0 1 2 3; 0 1 2 3]; 

  i43 = [0 0 0; 1 1 1; 2 2 2; 3 3 3];
  j43 = [0 1 2; 0 1 2; 0 1 2; 0 1 2]; 

       
end

function [jb1, jb2] = dcay_dq(q1, q2, cay, b)

% Jacobian entries - rotations wrt q
% Formulas from the bottom of page 9 in Rob's notes

  jb1 = zeros(3,4);
  jb2 = zeros(3,4);
  for i = 1:3
      B = reshape(b(i,:,:),4,4);
      jb1(i,:) = 1/(q2'*q1) * (-10*B - cay(i)*eye(4)) * q2;
      jb2(i,:) = 1/(q2'*q1) * ( 10*B - cay(i)*eye(4)) * q1;     
  end    
end

function [jb3, jb4] = dtrans_doq(q1, q2, o1, o2, b, iph)

% Jacobian entries - translations wrt o and q.
% Formulas from the top of page 10 in Rob's notes.
% Phosphate translations are expressed in the base (not middle) frame -
% iph indicates if we deal with a phosphate
  
  if iph 
      quat = q1;
  else
      quat = q1 + q2;
  end
  
  D = compute_ds(quat);
  jb3 = D';
  
  jb4 = zeros(3,4);
  DD = compute_dds(quat, D, b);
  
  for i = 1:3
    jb4(i,:) =  reshape(DD(i,:,:),3,4)' * (o2 - o1);   
  end
end

function [hb11, hb22, hb12] = ddcay_dqdq(q1, q2, cay, b)

% Hessian entries - rotations wrt qs
% Formulas from the page 12 in Rob's notes

  hb11 = zeros(3,4,4);
  hb22 = zeros(3,4,4);
  hb12 = zeros(3,4,4);
  for i = 1:3
      B = reshape(b(i,:,:),4,4);
      hb11(i,:,:) = (q2'*q1)^(-2) * (2*cay(i)*q2*(q2') + 10*(B*q2)*(q2') + 10*q2*(B*q2)');
      hb22(i,:,:) = (q2'*q1)^(-2) * (2*cay(i)*q1*(q1') - 10*(B*q1)*(q1') - 10*q1*(B*q1)');    
      hb12(i,:,:) = (q2'*q1)^(-2) * (2*cay(i)*q1*(q2') - 10*(B*q1)*(q2') + 10*q1*(B*q2)') ...
                                    + (q2'*q1)^(-1)*(10*B - cay(i)*eye(4));   
  end    
end

function [hbo, hbqq] = ddtrans_ddoq(q1, q2, o1, o2, b, iph)

% Hessian entries - translations wrt os and qs
% Formulas from the top of page 13 in Rob's notes
% Phosphate translations are expressed in the base (not middle) frame -
% iph indicates if we deal with a phosphate
  
  if iph 
      quat = q1;
  else
      quat = (q1 + q2); 
  end

  hbo = zeros(3,4,3);
  hbqq = zeros(3,4,4);
  
  dirs   = compute_ds(quat);
  ddirs  = compute_dds(quat, dirs, b);
  ddirs2 = compute_dds2(quat, dirs, ddirs, b);
  
  for j = 1:3
      
      hbo(j,:,:) = reshape(ddirs(j,:,:),3,4)';
      for k = 1:3
         hbqq(j,:,:) = reshape(hbqq(j,:,:),4,4)  ...
                        + (o2(k) - o1(k)) * reshape(ddirs2(j,k,:,:),4,4);                             
      end
  end    
end

function dirs = compute_ds(qvec)

% Compute the rotation matrix dirs corresponding to the quaternion qvec

  dirs = zeros(3,3);
  q1 = qvec(1); q2 = qvec(2); q3 = qvec(3); q4 = qvec(4);
  qnormsq = q1^2+q2^2+q3^2+q4^2;
  dirs(:,1) = (1/qnormsq)*[q1^2-q2^2-q3^2+q4^2; 2*q1*q2+2*q3*q4;      2*q1*q3-2*q2*q4];
  dirs(:,2) = (1/qnormsq)*[2*q1*q2-2*q3*q4;     -q1^2+q2^2-q3^2+q4^2; 2*q1*q4+2*q2*q3];
  dirs(:,3) = (1/qnormsq)*[2*q1*q3+2*q2*q4;     -2*q1*q4+2*q2*q3;     -q1^2-q2^2+q3^2+q4^2];
  
end

function ddirs = compute_dds(qvec, dirs, b)

% Calculate the first order partial derivatives of a 
% rotation matrix with respect to the quaternion components. 
% The inputs are a 4-element quaternion 'qvec' representing a rotation, 
% the rotation matric 'dirs' and the constant matrix 'b'.
% The output `ddirs` is a 3x3x4 tensor containing the partial derivatives.

ddirs = zeros(3,3,4);

for i = 1:3
    
  ii = [1 2 3 1 2];  
  i2 = ii(i+1);
  i3 = ii(i+2);  
  
  B2 = reshape(b(i2,:,:),4,4);
  B3 = reshape(b(i3,:,:),4,4); 
    
  ddirs(i,:,:) = dirs(:,i2)*(B3*qvec)' - dirs(:,i3)*(B2*qvec)';
    
end    
ddirs = 2*ddirs/(qvec'*qvec);

end

function ddirs2 = compute_dds2(qvec, dirs, ddirs, b)

% Calculate the second order partial derivatives of a 
% rotation matrix with respect to the quaternion components. 
% The inputs are a 4-element quaternion 'qvec', the constant matrix 'b',
% the rotation matrix 'dirs' and the first order partial derivatives of a 
% rotation matrix with respect to the quaternion 'ddirs'.
% The output 'ddirs2' is a 3x3x4x4 tensor containing the partial derivatives.

ddirs2 = zeros(3,3,4,4);

for i = 1:3
  
  ii = [1 2 3 1 2];  
  i2 = ii(i+1);
  i3 = ii(i+2);  
  
  B2 = reshape(b(i2,:,:),4,4);
  B3 = reshape(b(i3,:,:),4,4);
    
  for k = 1:3
      ddirs2(i,k,:,:) = dirs(k,i2)*B3 - dirs(k,i3)*B2 ... 
                        + (B3*qvec)*reshape(ddirs(i2,k,:),1,4) ...
                        - (B2*qvec)*reshape(ddirs(i3,k,:),1,4) ...
                        - reshape(ddirs(i,k,:),4,1)*qvec';                               
  end 
end  
ddirs2 = 2*ddirs2/norm(qvec)^2;   

end

function cay = cayley_from_q(q1, q2, b)

  q1 = q1/norm(q1);
  q2 = q2/norm(q2);
   
  cay = [0 0 0]';
  for i=1:3   
     Bi = reshape(b(i,:,:),4,4); 
     cay(i) = 10*q2'*Bi*q1/(q2'*q1);  
  end     

end
