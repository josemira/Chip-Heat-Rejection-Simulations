function steadyheat_equal_dxdy 
%-------------------------------------------------------------------------- 
% solving a steady-state conduction of a rectangular fin 
% Finite difference method with dx=dy 
%-------------------------------------------------------------------------- 
% physical dimension and parameters 
lx=0.012;  
ly=0.004; 
Tb=85; 
h=100; 
k=180; 
Tinf=20; 
%-------------------------------------------------------------------------- 
% numerical domain 
mx=48;  
my=16;  
nx=mx+1;  
ny=my+1; 
dx=lx/mx;  
dy=ly/my; 
if abs(dx-dy)>0.000001     
    fprintf('dx not equal to dy. Make sure you have modified the code\n'); 
end
x=linspace(0,lx,nx);
y=linspace(0,ly,ny); 
%-------------------------------------------------------------------------- 
% allocate u and initize it to be u_wall 
T=zeros(nx,ny);  
T(:,:)=Tb; 
%-------------------------------------------------------------------------- 
% const boundary conditions 
%left b.c. 
T(1,:)=Tb;   
%-------------------------------------------------------------------------- 
Tnew=T; 
%set(gca,'Tlim',[0 Tb]); 
%surf(x,y,T'); shading interp; 
%fprintf('Pressing any key to start...\n'); 
%pause; 
err_tol=0.0001; 
err=2*err_tol; 
iter=0;  
%-------------------------------------------------------------------------- 
% main loop 
while (err>err_tol) && iter<200000 
    %% i=1 is the left boundary, T=Tb, no need to update it     
    for i=2:mx         
        j=1;   % bottom boundary              
        Tnew(i,j)=(T(i-1,j)+ T(i+1,j)+2*T(i,j+1)+2*h*dx/k*Tinf)/(4+2*h*dx/k);  

        for j=2:my  % iternal nodes             
             Tnew(i,j)=(T(i-1,j)+T(i+1,j)+T(i,j-1)+T(i,j+1))/4;         
        end;

        j=my+1;  % top boundary         
        Tnew(i,j)=(T(i-1,j)+T(i+1,j)+T(i,j-1)*2+2*h*dx/k*Tinf)/(4+2*h*dx/k);          
    end;     

    i=mx+1; % right boundary         
        j=1; % bottom-right corner node                
        Tnew(i,j)=(T(i-1,j)+ T(i-1,j)+2*T(i,j+1)+2*h*dx/k*Tinf)/(4+2*h*dx/k);   

        for j=2:my  % righ boundary            
            Tnew(i,j)=(T(i-1,j)+T(i-1,j)+T(i,j-1)+T(i,j+1))/4;             
        end;         
        j=my+1; % top-right corner node         
        Tnew(i,j)=(T(i-1,j)+ T(i-1,j)+T(i,j-1)*2+2*h*dx/k*Tinf)/(4+2*h*dx/k); 
        
    err=sum(sum(abs(Tnew-T)));  
    
    T=Tnew;     
    iter=iter+1; 
end;  
%-------------------------------------------------------------------------- 
% output results 
fprintf('\nIter=%d, sum error=%f\n', iter,err); 
set(gca,'zlim',[0 Tb]); 
surf(x,y,T');shading interp;   

fileID=fopen('fin.txt','w'); 
for j=1:ny    %size(T,1)     
    fprintf(fileID, '%f\t', T(:,j));     
    fprintf(fileID,'\n'); 
end; 
fclose(fileID); 
%-------------------------------------------------------------------------- 