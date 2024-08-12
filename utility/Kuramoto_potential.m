function P = Kuramoto_potential(fai,param)
% ¼ÆËã Kuramoto ÊÆº¯Êý
fai = squeeze(fai);
A = param.Mat;
alpha = param.alpha;
N = param.node;

while ~(size(fai,1) == N)
    fai = fai';
end

if size(fai,2) == 1
    if isreal(param.Mat)
        P = -0.5 * alpha * sum(sum(A.*cos(fai - fai'))) / N; 
    else
        P = -0.5 * alpha * sum(sum(abs(A).*cos(fai - fai' - angle(A)))) / N;
    end
    
else
    Fai = repmat(fai - fai', 1, 1, size(fai,2));
    P = -0.5 * alpha * sum(sum(A.*cos(Fai))) / N;
end