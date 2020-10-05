function [hess,varargout]=hessian(f,p,varargin)

if nargout>2
    error('maximum # output arguments is 2');
end

fmt='%5.0f'; % printing format

nparam=length(p);
% Compute hessian element by element
ndx=6; % #purturbations 
h0=exp(-(0.05:0.01:0.1));
% h0=exp(-(0.01:0.01:0.06));

hess=zeros(nparam,nparam);
hessdiag=zeros(ndx,1);

dxscale=ones(nparam,1);
dxscale=sparse(1:nparam,1:nparam,dxscale,nparam,nparam);

fx=feval(f,p,varargin{:}); % evaluate objective function @ parameter values

disp(' ');
disp(' ');
disp('Computing Hessian...');
disp('---------------------------------------------');
disp(sprintf(['  diagonal elements    :   ' fmt ],nparam));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagonal elements of Hessian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nparam % for each diagonal element
    h=dxscale(:,i)*h0; 
    
    for j=1:ndx % for each purturbation, compute numerical differentiation
        
        % forward point
        dx=p+h(:,j);
        % backward point
        dy=p-h(:,j);
        % evaluate function at forward and backward point
        fdx=feval(f,dx,varargin{:});
        fdy=feval(f,dy,varargin{:});
        % Hessian at fdx and fdy value
        hessdiag(j)=-(2*fx-fdx-fdy)/(h(i,j)^2);
        
    end
    
    % diagonal value of Hessian is the average of the middle 2 purturbation
    % values 
    hess(i,i)=0.5*(hessdiag(3)+hessdiag(4));
    disp(sprintf(['                      ' fmt ],i));
    
end

disp(' ');
disp(sprintf(['  off-diagonal elements  :   ' fmt ],nparam*(nparam-1)/2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off-diagonal elements of Hessian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorij=[];
zz=1;

for i=1:nparam
    hi=dxscale(:,i)*h0;
    
    for j=(i+1):nparam
        hj=dxscale(:,j)*h0;
        
        for k=1:ndx
            
            % forward point
            dx=p+hi(:,k);
            % backward point
            dy=p-hj(:,k);
            % cross point
            dxdy=p+hi(:,k)-hj(:,k);
            % evaluate fucntion at points dx,dy and dxdy
            fdx=feval(f,dx,varargin{:});
            fdy=feval(f,dy,varargin{:});
            fdxdy=feval(f,dxdy,varargin{:});
            % Hessian at fdx,fdy,fdxdy values
            hessdiag(k)=-(fx-fdx-fdy+fdxdy)/(hi(i,k)*hj(j,k));
            
        end
        
        % off-diagonal (i,j) value of Hessian is the average of middle 2
        % purturbation values
        hess(i,j)=0.5*(hessdiag(3)+hessdiag(4));
        
        % calculating correlation among (i,j) values of Hessian
        if hess(i,i)==0 || hess(j,j)==0
            corrij=0;
        else
            corrij=hess(i,j)/sqrt(hess(i,i)*hess(j,j));
        end
        
        if abs(corrij) > 0.98 % correlation too big, change value
            hess(i,j)=0.9*sqrt(hess(i,i)*hess(j,j));
            errorij=[errorij; i j corrij];
        elseif abs(corrij) < 0.005 % correlation too small, make it 0
            hess(i,j)=0;
        end
        
        hess(j,i)=hess(i,j);
        
        if mod(zz,5)==0
            disp(sprintf(['                               ' fmt ],zz));
        end
        zz=zz+1;
        
    end
    
end

hess=real(hess);
varargout={errorij};


end
