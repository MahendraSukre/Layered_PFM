function PlotMatProp(X,T,Modulus)

% Variable Description:
%           X - The nodal coordinates of the mesh
%           -----> coordinates = [node X Y ] 
%           T - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]      
%           component - The components of stress whose profile to be plotted
%                     - 1 : Sigma_x, 2 : Sigma_y, 3 = Sigma_xy
%--------------------------------------------------------------------------
 nOfElements = size(T,1);
 Moduli = [Modulus Modulus Modulus Modulus];
 
 for iel = 1:nOfElements

%   get the current element connectivity
    econ = T(iel,:);
    
%   get the nodal coordiantes
    nds = X(econ,:);
    fill(nds(:,1), nds(:,2), Moduli(iel,:)); hold on
 end
 
colorbar