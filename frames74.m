function bp_level = frames74(shapes)  

    % Reconstruction of DNA base-pair, base and phosphate positions 
    % and orientations, starting from the middle, 74th, base-pair 
    % (nucleosome dyad point).
    %
    % The coordinate system is positioned in the center of the nucleosome,
    % with z axis going up vertically, x axis going from z towards the 
    % dyad point and y = cross(z,x).
    %
    % Input: shapes - a cgNA+ configuration vector for a 147bp DNA sequence,
    % the length of shapes is 3510.
    
    nbp = 147;
    
    % Absolute coordinates of the dyad point
    
    R0 = [-1 0 0; 0 0 1; 0 1 0];
    r0 = [40; 0; 0];
    
    G = R0;
    q = r0;
    
    % Relative coordinates of the oligomer
    [eta, w, etapW, wpW, u, v, etapC, wpC] = vector2shapes(shapes);
    
    bp_level = InitializeStruct(nbp) ;
    
    for i = 74:nbp % second half of DNA (forward reconstruction)
           
        bp_level(i).R = G ;
        bp_level(i).r = q ;
        
        % base pair:
        r = cay(eta(i,:));
        Gw = G * w(i,:)';  
 
        % complementary strand 
        bp_level(i).Rc = G * (sqrtm(r))'; 
        bp_level(i).rc = q - 0.5 * Gw;
        
        % main strand
        bp_level(i).Rw = bp_level(i).Rc * r;
        bp_level(i).rw = bp_level(i).rc + Gw;

        if i < nbp
            ru = cay(u(i,:));
            sqrtru = sqrtm(ru);
            H = G * sqrtru;
            % next base pair:
            G = G * ru;
            q = q + H * v(i,:)';   
        end  
    end
    
    G = R0;
    q = r0;
       
    for i = 74:-1:1 % first half of DNA (backward reconstruction)
            
        bp_level(i).R = G ;
        bp_level(i).r = q ;
        
        % base pair:
        r = cay(eta(i,:));
        Gw = G * w(i,:)';  
 
        % complementary strand 
        bp_level(i).Rc = G * (sqrtm(r))'; 
        bp_level(i).rc = q - 0.5 * Gw;
        
        % main strand
        bp_level(i).Rw = bp_level(i).Rc * r;
        bp_level(i).rw = bp_level(i).rc + Gw;
          
        if i > 1
            ru = cay(-u(i-1,:));
            sqrtru = sqrtm(ru);
            H = G * sqrtru;
            % next base pair:
            G = G * ru;
            q = q - H * v(i-1,:)';   
        end  
    end
    
    for i = 1:nbp - 1
    
      bp_level(i+1).Rpw = bp_level(i+1).Rw*cay(etapW(i,:)) ;
      bp_level(i+1).rpw  = bp_level(i+1).rw + bp_level(i+1).Rw*wpW(i,:)' ;
      
      Rc = bp_level(i).Rc*diag([1,-1,-1]) ;
      
      bp_level(i).Rpc = Rc*cay(etapC(i,:)) ;
      bp_level(i).rpc  = bp_level(i).rc + Rc*wpC(i,:)' ;
      
    end
    
end

function bp_level = InitializeStruct(nbp)

  bp_level  = struct('R', [],'r',[], ...
             'Rw', [], 'rw', [], ...
             'Rc', [], 'rc', [], ...
             'Rpw', [], 'rpw', [], ...
             'Rpc', [], 'rpc', cell(nbp,1)) ;
end

function Q = cay(k)

    I = eye(3) ;
    alpha = 1/10 ;
    k = alpha*k ;
    X = [   0   -k(3)  k(2) ;
           k(3)   0   -k(1) ;
          -k(2)  k(1)   0 ] ;
    Q = (I+X)/(I-X) ;

end