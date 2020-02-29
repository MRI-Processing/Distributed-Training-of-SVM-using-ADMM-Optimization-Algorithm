rand('seed', 0);
randn('seed', 0);

n = 2;
m = 200;
N = m/2;
M = m/2;

rho=1;
alpha=1;

% positive examples
Y = [1.5+0.9*randn(1,0.6*N), 1.5+0.7*randn(1,0.4*N);
2*(randn(1,0.6*N)+1), 2*(randn(1,0.4*N)-1)];

% negative examples
X = [-1.5+0.9*randn(1,0.6*M),  -1.5+0.7*randn(1,0.4*M);
2*(randn(1,0.6*M)-1), 2*(randn(1,0.4*M)+1)];

x = [X Y];                           %---------->total dataset
y = [ones(1,N) -ones(1,M)];          %----------> class labels
A = [ -((ones(n,1)*y).*x)' -y'];   
xdat = x';
xdat(1:10,:)
lambda = 1.0;


% p = zeros(1,m);
% p(y == 1)  = sort(randi([1 10], sum(y==1),1));
% p(y == -1) = sort(randi([11 20], sum(y==-1),1));
% 
% 
% [m, n] = size(A);
% N = max(p);
% 
% for i = 1:N
%     tmp{i} = A(p==i,:);
% end
% A = tmp;
% 
% x = zeros(n,N);
% z = zeros(n,N);
% u = zeros(n,N);
% 
% for k = 1:2
% 
% 	% x-update
%     for i = 1:N
%         cvx_begin quiet
%             variable x_var(n)
%             
%             minimize ( sum(pos(A{i}*x_var + 1)) + rho/2*sum_square(x_var - z(:,i) + u(:,i)) )
%         cvx_end
%         x(:,i) = x_var;  
%     end
%     xave = mean(x,2);
% 
%     % z-update with relaxation
%     zold = z;
%     x_hat = alpha*x +(1-alpha)*zold;
%     z = N*rho/(1/lambda + N*rho)*mean( x_hat + u, 2 );
%     z = z*ones(1,N);
% 
%     % u-update
%     u = u + (x_hat - z);
%  
%  end
% 

