function zvec = z_from_w_bp_28ph_o_center(wvec)

% Computes a vector z of hybrid coordinates from the cgNA+ coordinate
% vector wvec

global indc
global indw
global nn
global nc

bp_level  = frames74(wvec);
 
%indc =  [6 16 25 36 46 56 66 76 87 97 108 118 127 137];
%indw =  [11 21 30 40 51 61 72 82 92 102 112 123 132 142];

nbp = (length(wvec)+18)/24;

qbc = zeros(nbp,4);
obc = zeros(nbp,3);
qbw = zeros(nbp,4);
obw = zeros(nbp,3);
opc = zeros(nbp,3);
opw = zeros(nbp,3);
qbp = zeros(nbp,4);
obp = zeros(nbp,3);
   

for i=1:nbp
    
    if ((i>1)&&(q1'*rot_to_quaternion(bp_level(i).R)<0)) || ((i==1)&&(rot_to_quaternion(bp_level(i).R)'*rot_to_quaternion(bp_level(74).R)<0))
       q1 = - rot_to_quaternion(bp_level(i).R);
    else
       q1 = rot_to_quaternion(bp_level(i).R);
    end    
    
    qbp(i,:)  = q1;
    obp(i,:)  = bp_level(i).r;
    
    if sum(i == indc)
 
      qbc(i,:)  = rot_to_quaternion(bp_level(i).Rc);  
      obc(i,:)  = bp_level(i).rc; 
      opc(i,:)  = bp_level(i).rpc;
      
      if (qbp(i,:)*qbc(i,:)'<0)
         qbc(i,:) = -qbc(i,:);
      end    
      
    elseif sum(i == indw)
      
      qbw(i,:)  = rot_to_quaternion(bp_level(i).Rw);
      obw(i,:)  = bp_level(i).rw;
      opw(i,:)  = bp_level(i).rpw;
      
      if (qbp(i,:)*qbw(i,:)'<0)
         qbw(i,:) = -qbw(i,:);
      end   
      
    end     
end

nc = length([indc indw]);
zvec = zeros(nbp*25-12+nc,1);

zvec(1:3) = obp(1,:);  
zvec(4:7) = qbp(1,:);
zvec(8:19)  = wvec(1:12);
k2 = 19;

for i = 2:nbp
   
   zvec(k2+1:k2+3) = obp(i,:);  
   zvec(k2+4:k2+7) = qbp(i,:);
   k1 = (i-1)*24; 

   if i<nbp
      if sum(i == indc) 
         zvec(k2+8:k2+10)  = obc(i,:); 
         zvec(k2+11:k2+14) = qbc(i,:);
         zvec(k2+15:k2+17) = opc(i,:); 
         zvec(k2+18:k2+20) = wvec(k1+7:k1+9); % C phosphate rotations
         zvec(k2+21:k2+26) = wvec(k1-5:k1);   % W phosphate rotations & translations
         k2 = k2+26;
      elseif sum(i == indw)
         zvec(k2+8:k2+10)  = obw(i,:); 
         zvec(k2+11:k2+14) = qbw(i,:);
         zvec(k2+15:k2+17) = opw(i,:); 
         zvec(k2+18:k2+20) = wvec(k1-5:k1-3);  % W phosphate rotations
         zvec(k2+21:k2+26) = wvec(k1+7:k1+12); % C phosphate rotations & translations
         k2 = k2+26;
      else 
         zvec(k2+8:k2+19)  = wvec(k1+1:k1+12);
         zvec(k2+20:k2+25) = wvec(k1-5:k1);
         k2 = k2+25;   
      end
   else
      zvec(k2+8:k2+13)  = wvec(k1+1:k1+6); 
      zvec(k2+14:k2+19) = wvec(k1+19-24:k1+24-24); 
   end    
end  

% bp number 74 (the dyad point) is fixed
zvec = zvec([1:nn,nn+8:end]);
  
end

