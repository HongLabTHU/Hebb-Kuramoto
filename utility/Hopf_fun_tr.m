function xdot = Hopf_fun_tr(t ,x, param)
Noise = param.Noise;
M = param.Mat;
node = param.node;
alpha = param.alpha;
dt = param.dt;

a = param.a;
b = param.b;

% sigma = param.sigma;

% node 对应
x0 = x(1:node);
x1 = x(node+1 : 2*node);

%
% x1(x1>pi) = x1(x1>pi) - 2*pi;
% x1(x1<-pi) = x1(x1<-pi) + 2*pi;

dig = zeros(1, node)';  % linspace(-1,1,6)';
x1m = M*exp(1i*x1) / node;
inp = inputfunc(t, param) + alpha * abs(x1m) .* sin(angle(x1m)-x1-dig);

if ~isempty(Noise)
    inp = inp + Noise(:, floor(t/dt)+1);
end

dx0 = a*x0 - x0.^3 ;
dx1 = b + inp ;  % + normrnd(0,sigma,node,1);  % kuramoto形式

xdot=[dx0; dx1];


end

function out = Siga(v)
% 激活函数
vmax = 5;  % Hz  Maximum firing rate
v0 = 6;  % mV Potential at half of maximum firing rate
r = 0.56;  %  mV  Slope of sigmoid function at v0 

out = vmax ./ (1 + exp(r*(v0-v)));
end


function out = inputfunc(t, param)
% 输入函数
inputtime = param.inputtime;    % 输入时间
inputpulse = param.inputpulse;     % 输入大小
inputstd = param.inputstd;
out = inputpulse * exp(-(t-inputtime)^2/(inputstd^2));
end