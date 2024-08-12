function xdot = Hopf_fun(t ,x, param)
% Noise = param.Noise;
dt = param.dt;
M = param.Mat;
node = param.node;
alpha = param.alpha;

a = param.a;
b = param.b;

% node ∂‘”¶
x0 = x(1:node);
x1 = x(node+1 : 2*node);


dx0 = 0*x0;  %a.*x0 - x0.^3 ;
dx1 = b + alpha * diag(M * sin(x1 - x1'))/node;

xdot=[dx0; dx1];


end
