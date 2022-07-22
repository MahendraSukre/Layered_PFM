function [Strain, Stress, Modulus] = Strain_Stress(X,T,Xref,Tref,refinedElements,E_elems,nu_elems,ux,uy,MatTyp,We,m, M)

nOfElements = size(T,1);
standardElements = setdiff(1:nOfElements,refinedElements);
nOfNodes = size(X,1) + size(Xref,1);
Stress = zeros(nOfElements,size(T,2),3);
Strain = zeros(nOfElements,size(T,2),3);
Modulus = zeros(nOfElements,1);

% Volume integrals (loop in elements)
for i = 1:length(standardElements)
    iElem = standardElements(i);
    Te = T(iElem,:) ;
    Xe = X(Te,:);
    switch MatTyp
   case 1
        E = E_elems(iElem);
   case 2
      if Xe(1,2) < 0
        E = E_elems(iElem);
      else
        E = We*E_elems(iElem);
      end
   case 4
       X1 = Xe(:,1); Y1 = Xe(:,2);
       polyin = polyshape(X1' , Y1', 'Simplify', false);
       [xc,yc] = centroid(polyin);
       yb = m*xc;
     if yc < yb
        E = E_elems(iElem);
     else
        E = We*E_elems(iElem);
     end
   otherwise
      error('Material type not defined')
    end
    
    Modulus(iElem,1) = E;
    nu = nu_elems(iElem);
    U = [ux(Te) ; uy(Te)];
    [strain, stress] = computeVolumeMatrices(Xe,E,nu,U);
    Stress(iElem,:,:) = stress;
    Strain(iElem,:,:) = strain;
end

% About refined elements

for i = 1:length(refinedElements)
    iElem = refinedElements(i);
    Te = Tref(i,:);
    Tc = Te([1,M,M^2,M*(M-1)+1]);
    Xe = Xref(Te,:); 
    Xc = Xe([1,M,M^2,M*(M-1)+1],:);
    switch MatTyp
   case 1
        E = E_elems(iElem);
   case 2
    if Xe(2,1) <= 0
        E = E_elems(iElem);
    else
        E = We*E_elems(iElem);
    end  
   case 4
       X1 = Xe(:,1); Y1 = Xe(:,2);
       polyin = polyshape(X1' , Y1', 'Simplify', false);
       [xc,yc] = centroid(polyin);
       yb = m*xc;
     if yc < yb
        E = E_elems(iElem);
     else
        E = We*E_elems(iElem);
     end
    otherwise
      error('Material type not defined')
    end
    
    Modulus(iElem,1) = E;
    nu = nu_elems(iElem);
    ind = [Tc,Tc + nOfNodes];
    U = [ux(Tc) ; uy(Tc)];
    [strain, stress] = computeVolumeMatrices(Xc,E,nu,U);
    Stress(iElem,:,:) = stress;
    Strain(iElem,:,:) = strain;
end  
end

function [strain, stress] = computeVolumeMatrices(Xe,E,nu,U)
        
        strain =  zeros(4,3); stress = zeros(4,3);
        
        % Nodal position in master element
        vertex = [0 0; 1 0; 1 1; 0 1];   
        nds = Xe;
        D = E/((1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 0.5*(1-2*nu)];

for nn = 1:length(vertex)
        
        % get the current integration point
        pt = vertex(nn,:);
        
        % get the shape functions and its derivatives
        [~,dNdxi] = lagrange_basis('Q4',pt,2)  ;
       
        % get the jacobian
        J = (dNdxi'*nds);
        
        % get the derivatives in the global space
        dNdx = (J\dNdxi')';
        
        Be  = [dNdx(1,1) dNdx(2,1) dNdx(3,1) dNdx(4,1) 0 0 0 0;
                0 0 0 0 dNdx(1,2) dNdx(2,2) dNdx(3,2) dNdx(4,2);
                dNdx(1,2) dNdx(2,2) dNdx(3,2) dNdx(4,2) dNdx(1,1) dNdx(2,1) dNdx(3,1) dNdx(4,1)];
     
        strain0 = Be*U;
        strain(nn,:) = strain0;
        stress(nn,:) = D*strain0;
end 
end