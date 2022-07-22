function u_D = dirichletFunction(X,k,test)

 x = X(:,1); 
 y = X(:,2);
 
u_D = zeros(length(x),2); %[ux,uy]
if (test==1)
    dirichletA = abs(y-0.5)<1.e-6;
    u_D(dirichletA,1) = k;    
elseif (test==5)
    u_D(y > 80,2) = k;   
elseif (test == 6)
    dirichletA = abs(y - 1)==0;  % top face
    u_D(dirichletA,2) = k;             % Disp. BC
elseif (test == 7)
    dirichletA = abs(x - 50)==0;  % top face
    u_D(dirichletA,1) = k ;
elseif (test == 8)
    dirichletA = abs(x - 20)==0;  % top face
    u_D(dirichletA,1) = k;
end


