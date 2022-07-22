function uCCD = dirichletCondition(X,Xref,nodesCCDStd,nodesCCDRef,k,test,BC)

             
XdirichletStd = X(nodesCCDStd,:);
XdirichletRef = Xref(nodesCCDRef,:);
fix = min(X(:,1));

uCCDStd = dirichletFunction(XdirichletStd,k,test);

if(nodesCCDRef)
    
    uCCDRef = dirichletFunction(XdirichletRef,k,test);
    
    if BC == 1         % 1: Right edge free in Y, 
       BC_Ref = Xref(nodesCCDRef,:)== fix;
       uCCD = [uCCDStd(:,1); uCCDRef(:,1); uCCDStd(1:size(uCCDStd,1)/2,2); uCCDRef(BC_Ref,2)]; %[uxS;uxR;uyS;uyR]
    
    elseif BC == 2     % 2: L&R both edges free in Y except at (-20,0), else fixed in Y
       BC_Std = X(nodesCCDStd,1) == fix & X(nodesCCDStd,2) == 0; 
       BC_Ref = Xref(nodesCCDRef,1) == fix & Xref(nodesCCDRef,2) == 0;
       uCCD = [uCCDStd(:,1); uCCDRef(:,1); uCCDStd(BC_Std,2); uCCDRef(BC_Ref,2)];
   
    elseif BC == 3     % three point bending
       BC_Std = abs(X(nodesCCDStd,1))== 10 & X(nodesCCDStd,2)== -3;
       BC_Ref = abs(Xref(nodesCCDRef,1))== 10 & Xref(nodesCCDRef,2)== -3;
       uCCD = [uCCDStd(BC_Std,1); uCCDRef(BC_Ref(1),1); uCCDStd(BC_Std,2); uCCDRef(BC_Ref,2)];
    
    else 
       uCCD = [uCCDStd(:,1); uCCDRef(:,1); uCCDStd(:,2); uCCDRef(:,2)];
    end
    
else
    if BC == 1              % 1: Right edge free in Y, 
       uCCD = [uCCDStd(:,1); uCCDStd(1:size(uCCDStd,1)/2,2)];
   
    elseif BC == 2          % 2: L&R both edges free in Y except at (-20,0), else fixed in Y
       BC_Std = X(nodesCCDStd,1) == fix & X(nodesCCDStd,2) == 0;
       uCCD = [uCCDStd(:,1); uCCDStd(BC_Std,2)];
       
    elseif BC == 3          % three point bending
       BC_Std = find(abs(X(nodesCCDStd,1))== 10 & X(nodesCCDStd,2)== -3);
       uCCD = [uCCDStd(BC_Std(1),1); uCCDStd(BC_Std,2)];
    
    else
       uCCD = [uCCDStd(:,1); uCCDStd(:,2)];
    end
    
end
end