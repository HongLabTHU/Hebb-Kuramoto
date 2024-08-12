function [res, Aff] = wocc_projection(v, tem, mode)
% 计算wocc投影
% v: 【time, node】
% tem: 【node, 4】
% mode=1：
% mode=2：归一化

N_v = sqrt(sum(abs(v).^2, 2));
N_t = sqrt(sum(abs(tem).^2, 1));
Aff = v * conj(tem) ./ (N_v * N_t);

% if mode == 2
%     Aff = v * conj(tem) ./ (abs(v) * abs(tem));
% end

res = zeros(size(v,1), 2);
res(: , 1) = abs(Aff(:, 1)) - abs(Aff(:, 2));
res(: , 2) = abs(Aff(:, 3)) - abs(Aff(:, 4));

if mode == 2
    res = zeros(size(v,1), 2);
    res(: , 1) = abs(Aff(:, 1)).^2 - abs(Aff(:, 2)).^2;
    res(: , 2) = abs(Aff(:, 3)).^2 - abs(Aff(:, 4)).^2;
end

end