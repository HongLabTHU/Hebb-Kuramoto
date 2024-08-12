function J = Kuramoto_Jakobi(fai,param, mode)
% 注意该函数只适用于计算不变流形的雅可比矩阵
M = param.Mat;
alpha = param.alpha;
N = param.node;
% 
% D = abs(A*exp(1i*fai));
% 
% J = alpha * ( A.*(cos(fai')./cos(fai)) - diag(D) );
if mode == 1
vars = [];
for i = 1:N  % length(fai)
    vars = [vars ' s' num2str(i) ''];
end
eval(['syms ', vars]);

eval(['x1 = [' vars '];']);
x1 = x1';
f = alpha * diag(M * sin(x1 - x1'))/N;
J = eval(['jacobian(f, [', vars, '])']);

for i = 1:N
    eval(['s' num2str(i) ' = ' 'fai(' num2str(i) ');'])
end

J = eval(J);
end

if mode == 2
    Q = M .* cos(fai - fai');
    J = (alpha/N) * (Q - diag(sum(Q,2)));
end

end
