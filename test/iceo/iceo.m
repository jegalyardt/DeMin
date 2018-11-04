% function [f,low,upp,fmin,xmin]=iceo(prob,x)
% test functions from the 
% First International Contest on Evolutionary Optimization
%
% prob:		problem number 
% x:		argument
% 
% f:		function value (when nargout==1)
% 		value to reach (when nargout>1)
% low: 		vector of lower bounds of problem definition
% upp:		vector of upper bounds of problem definition
% fmin:		global minimum (or [])
% xmin:		global minimizer (or [])
% 
function [f,low,upp]=iceo(prob,x)

if prob==1,
  % The  sphere model
  % separable, quadratic, unimodal
  if nargout==1,
    f=(x-1)'*(x-1);
  else
    f=1e-6;
    fmin=0;
    xmin=1+0*x;
    low=-5*ones(size(x));
    upp=5*ones(size(x));
  end;
  return;
end;

if prob==2,
  % Griewank's function
  % separable at minimizer of partial objective
  if nargout==1,
    d=4000;
    prod=1;
    for i=1:size(x,1), prod=prod*cos((x(i)-100)/sqrt(i)); end;
    f=(x-100)'*(x-100)/d-prod+1;
  else
    f=1e-4;
    fmin=0;
    xmin=100+0*x;
    low=-600*ones(size(x));
    upp=600*ones(size(x));
  end;
  return;
end;

if prob==3,
  % Shekel's foxholes
  % rational, nonseparable
  [A,c]=iceo3d;
  if nargout==1,
    dim=size(x,1);m=size(c,1);
    sqn=zeros(m,1);
    for i=1:m, 
      sqn(i)=(x-A(1:dim,i))'*(x-A(1:dim,i)); 
    end;
    f=-sum(1./(sqn+c));
  else
    if size(x,1)>10, error('dimension must be at most 10!'); end;
    f=-9;
    low=0*ones(size(x));
    upp=10*ones(size(x));
  end;
  return;
end;

if prob==4,
  % Michalewicz's function
  % separable multimodal
  if nargout==1,
    m=10;zweim=2*m;
    for i=1:size(x,1), fac(i)=sin(i*x(i)^2/pi); end;
    f=-sin(x)'*(fac'.^zweim);    
  else
    if size(x,1)==5, f=-4.687;
    elseif size(x,1)==10, f=-9.66; 
    else f=-0.966*size(x,1);
    end;
    low=0*ones(size(x));
    upp=pi*ones(size(x));
  end;
  return;
end;

if prob==5,
  % Langerman's function
  % exponential, nonseparable
  [A,c]=iceo5d;
  if nargout==1,
    dim=size(x,1);m=size(c,1);
    sqn=zeros(m,1);
    for i=1:m, sqn(i)=(x-A(1:dim,i))'*(x-A(1:dim,i)); end;
    f=-c'*(exp(-sqn/pi).*cos(sqn*pi));
  else
    if size(x,1)>10, error('dimension must be at most 10!'); end;
    f=-1.4;
    low=0*ones(size(x));
    upp=10*ones(size(x));
  end;
  return;
end;

% problems >5 not in original iceo test set
if prob=='rb',
  % Rosenbrock function
  % nonseparable
  if nargout==1,
    z1=x(2)-x(1)^2;
    z2=1-x(1);
    f=100*z1^2+z2^2;
  else
    if size(x,1)~=2, error('dimension must be 2!'); end;
    f=1e-4;
    low=-1000*rand*ones(size(x));
    upp=200*rand*ones(size(x));
  end;
  return;
end;

if prob=='mcubic',
  % 1D cubic modified by an oscillating term
  if nargout==1,
    f=(x+3).*x^2+3+3*cos(10*x); 
  else
    if size(x,1)~=1, error('dimension must be 1!'); end;
    f=0.2635;
    fmin=0.26093049893318;
    xmin=-0.30893243005125;
    low=-3.28664+2*rand;
    upp=3.28664-2*rand; 
  end;        
  return;
end;

if prob==' ',
  % 
  if nargout==1,
    f=0; 
  else
    if size(x,1)~=0, error('dimension must be 0!'); end;
    f=0;
    fmin=0;
    xmin=0;
    low=0;
    upp=0; 
  end;        
  return;
end;

error('function not defined')
  

