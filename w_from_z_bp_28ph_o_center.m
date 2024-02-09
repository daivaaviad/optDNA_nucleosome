function [wvec,zel] = w_from_z_bp_28ph_o_center(zvec)

% Computes the cgNA+ coordinate vector wvec from the vector 
% z of hybrid coordinates 

global indc
global indw
global nn
global nc

%indc =  [6 16 25 36 46 56 66 76 87 97 108 118 127 137];
%indw =  [11 21 30 40 51 61 72 82 92 102 112 123 132 142];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zlen  = length(zvec)+7;  

% bp number 74 (the dyad point) is fixed
R0 = [-1 0 0; 0 0 1; 0 1 0];
r0 = [40; 0; 0];
zvec = [zvec(1:nn); r0; rot_to_quaternion(R0); zvec(nn+1:end)];

nbp = (zlen+12-nc)/25;

qbp = zeros(nbp,4);
obp = zeros(nbp,3);
qbc = zeros(nbp,4);
obc = zeros(nbp,3);
qbw = zeros(nbp,4);
obw = zeros(nbp,3);
qpw = zeros(nbp,4);
opw = zeros(nbp,3);
qpc = zeros(nbp,4);
opc = zeros(nbp,3);

% The constant matrices Bi
b  = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3)  = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1)  = -1; b(3,3,4) = 1; b(3,4,3) = -1;  

B1 = reshape(b(1,:,:),4,4); 

obp(1,:) = zvec(1:3);
qbp(1,:) = zvec(4:7);
k2 = 19;

for i = 2:nbp
   obp(i,:) = zvec(k2+1:k2+3);
   qbp(i,:) = zvec(k2+4:k2+7);
   if sum(i == indc) 
      obc(i,:) = zvec(k2+8:k2+10);
      qbc(i,:) = zvec(k2+11:k2+14);
      opc(i,:) = zvec(k2+15:k2+17);
      %qpc(i,:) = zvec(k2+18:k2+21);
      k2 = k2 + 26;
   elseif sum(i == indw) 
      obw(i,:) = zvec(k2+8:k2+10); 
      qbw(i,:) = zvec(k2+11:k2+14);
      opw(i,:) = zvec(k2+15:k2+17);
      %qpw(i,:) = zvec(k2+18:k2+21); 
      k2 = k2 + 26;
   else 
      k2 = k2 + 25;
   end    
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inters_all = zeros(6,nbp-1);

for i=1:nbp-1
        
    R  = quaternion_to_rot(qbp(i,:) + qbp(i+1,:)); 
    inters_all(1:3,i) = cayley_from_q(qbp(i,:)', qbp(i+1,:)', b);
    inters_all(4:6,i) = R'*(obp(i+1,:)'- obp(i,:)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Intras and phosphates %%%%%%%%%%%%%%%%%

intras_all = zeros(6,nbp);
phosC = zeros(6,nbp);
phosW = zeros(6,nbp);

intras_all(1:6,1) = zvec(8:13); 
phosC(1:6,1) = zvec(14:19); 
k2 = 19;

for i=2:nbp
    
    if sum(i == indc) 
              
      R  = quaternion_to_rot(qbp(i,:)); 
      cay = cayley_from_q(qbc(i,:)', qbp(i,:)', b);
      intras_all(1:3,i) = cay*2/(1-norm(cay)^2/100);
      intras_all(4:6,i) = 2*R'*(obp(i,:) - obc(i,:))';
      
      qc1 = B1*qbc(i,:)';
      R  = quaternion_to_rot(qc1);
      phosC(4:6,i) = R'*(opc(i,:) - obc(i,:))';
      
      phosC(1:3,i) = zvec(k2+18:k2+20);
      phosW(1:6,i) = zvec(k2+21:k2+26); 
      
      k2 = k2 + 26;
      
    elseif sum(i == indw) 
  
      R  = quaternion_to_rot(qbp(i,:)); 
      cay = cayley_from_q(qbp(i,:)', qbw(i,:)', b);
      intras_all(1:3,i) = cay*2/(1-norm(cay)^2/100);
      intras_all(4:6,i) = 2*R'*(obw(i,:) - obp(i,:))';
      
      R  = quaternion_to_rot(qbw(i,:));
      phosW(4:6,i) = R'*(opw(i,:)'- obw(i,:)');
      
      phosW(1:3,i) = zvec(k2+18:k2+20);
      phosC(1:6,i) = zvec(k2+21:k2+26); 
        
      k2 = k2 + 26;
      
    elseif i<nbp 
        
      intras_all(1:6,i) = zvec(k2+8:k2+13); 
      phosC(1:6,i) = zvec(k2+14:k2+19); 
      phosW(1:6,i) = zvec(k2+20:k2+25); 

      k2 = k2 + 25;
    else  
      intras_all(1:6,i) = zvec(k2+8:k2+13); 
      phosW(1:6,i) = zvec(k2+14:k2+19); 
    end  

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Assembling the wvec   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvec = zeros(24*nbp-18,1);

wvec(1:18) = [intras_all(:,1); phosC(:,1); inters_all(:,1)];
wvec(end-11:end) = [phosW(:,nbp); intras_all(:,nbp)];

for i=2:nbp-1
    wvec(18+24*(i-2)+1:18+24*(i-2)+24) = [phosW(:,i); intras_all(:,i); phosC(:,i); inters_all(:,i)];
end

zel.qbp = qbp;
zel.obp = obp;
zel.qbc = qbc;
zel.obc = obc;
zel.qbw = qbw;
zel.obw = obw;
zel.qpc = qpc;
zel.opc = opc;
zel.qpw = qpw;
zel.opw = opw;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cay = cayley_from_q(q1, q2, b)

  q1 = q1/norm(q1);
  q2 = q2/norm(q2);
   
  cay = [0 0 0]';
  for i=1:3   
     Bi = reshape(b(i,:,:),4,4); 
     cay(i) = 10*q2'*Bi*q1/(q2'*q1);  
  end     

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






