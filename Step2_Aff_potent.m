%% memory two different pahse patterns
clear;clc

node = 36;   % 节点数
param.alpha = 200;
param.node = node;

F = Traveling_Mode;
grids = F.set_grids(6);
tem1 = F.rotational(grids,[0 0], 1, 1);  % 旋转型
tem2 = F.translational(grids,pi/2, 2);  % 前后型


fai = conj(tem1');
%^ A = cos(angle(fai)-angle(fai)');
A = [real(fai) imag(fai)] * pinv([real(fai) imag(fai)]);
param.Mat = A;

% 验证两个条件
% real(fai)'*imag(fai)
% norm(real(fai))
% norm(imag(fai))

% scatter(real(fai),imag(fai));
% hold on
% scatter(real(A*fai),imag(A*fai));

P = [];
Q = [];
U = [];
Fai = 0.5*pi*rand(200,param.node) + repmat(angle(conj(fai')),200,1);
Fai = exp(1i*Fai);
Fai = [conj(fai'); Fai];
for i = 1:201
    P(i) = Kuramoto_potential(angle(Fai(i,:)),param);
    U(i,1) = Fai(i,:)*fai;
    U(i,2) = Fai(i,:)*conj(fai);
    Q(i) = abs(Wfunc(Fai(i,:),fai.'))^2 + abs(Wfunc(Fai(i,:),fai'))^2;
    Q(i) = param.alpha * Q(i) / -2;
end

scatter(P(1), Q(1),10,[0 0.4470 0.7410], 'filled');hold on
plot(P,Q, 'color', [0.9290 0.6940 0.1250],'linewidth', 1.4);
scatter(P(2:200+1), Q(2:200+1),8,[0.4940 0.1840 0.5560], 'filled');  hold off
grid on
ylabel('-H($e^{i\theta}$, $e^{i\varphi_{R1}}$)', 'FontSize',14,'Interpreter', 'latex');
xlabel('P($\theta$, A)', 'FontSize',14,'Interpreter', 'latex');
title('memorize single pattern', 'FontSize',14,'Interpreter', 'latex');

%%
Rtem = [real(fai) imag(fai)];
Rtem = Rtem'*Rtem;
sl = 0.5:1:(size(Rtem, 1)-0.5);
subplot(1,2,1)
colormap cool;
imagesc(sl, sl, Rtem);
title('Orthogonal matrix')
grid on
% 在每个像素上添加数值
for i = 1:size(Rtem, 1)
    for j = 1:size(Rtem, 2)
        text(sl(j), sl(i), num2str(Rtem(i, j),4), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
colorbar

subplot(1,2,2)
Rtem = [real(fai) imag(fai)];
bar(sum(Rtem.^2).^0.5);
title('norm of  Re\phi  Im\phi ')


%% 验证正交性和范数条件
clear;clc

node = 36;   % 节点数
param.alpha = 200;
param.node = node;

F = Traveling_Mode;
grids = F.set_grids(6);
tem1 = F.rotational(grids,[0 0], 1, 1);  % 旋转型
tem2 = F.translational(grids,pi/2,2);  % 前后型

tem = [tem1; conj(tem1); tem2; conj(tem2)];
tem = conj(tem');

Rtem = [real(tem(:,1)) imag(tem(:,1)) real(tem(:,3)) imag(tem(:,3))];
Rtem = Rtem'*Rtem;
sl = 0.5:1:(size(Rtem, 1)-0.5);
subplot(1,2,1)
colormap cool;
imagesc(sl, sl, Rtem);
set(gca,'XTick',sl);
set(gca,'XTicklabel', {'Re(R1)','Im(R1)','Re(D1)','Im(D1)'});
set(gca,'YTick',sl);
set(gca,'YTicklabel', {'Re(R1)','Im(R1)','Re(D1)','Im(D1)'});
title('Orthogonal matrix')


% 在每个像素上添加数值
for i = 1:size(Rtem, 1)
    for j = 1:size(Rtem, 2)
        text(sl(j), sl(i), num2str(Rtem(i, j),4), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
colorbar

subplot(1,2,2)
Rtem = [real(tem(:,1)) imag(tem(:,1)) real(tem(:,2)) imag(tem(:,2))...
    real(tem(:,3)) imag(tem(:,3)) real(tem(:,4)) imag(tem(:, 4))];
bar(sum(Rtem.^2).^0.5);
title('norm of R1-R2-D1-D2')


%% 验证伪逆法 affinities 和 能量函数之间的关系
fai = tem(:,1);  % 基于模板4构造
fai2 = tem(:,3);
V1 = [real(fai) imag(fai)];
V2 = [real(fai2) imag(fai2)];
Mat = [0.5*V1 0.5*V2] * pinv([V1 V2]);  % 连接矩阵

param.Mat = Mat;

P = [];
Q = [];
U = [];
num = 1681;
Fai = [];

for j = 1:5
    Fai = [Fai; 0.2*j*pi*rand(200,param.node) + repmat(angle(fai.'),200,1)];
    Fai = [Fai; 0.2*j*pi*rand(200,param.node) + repmat(angle(fai2.'),200,1)];
end
Fai = exp(1i*Fai);
Fai = [conj(fai');conj(fai2'); Fai];

lamda = 0.48;  %0.01:0.02:0.99;
R2R = zeros(1,length(lamda));
for gg = 1:length(lamda)
    for i = 1:num+1
        P(i) = Kuramoto_potential(angle(Fai(i,:)),param);
        U(i,1) = abs(Wfunc(Fai(i,:),fai.'))^2 + abs(Wfunc(Fai(i,:),fai'))^2;
        U(i,2) = abs(Wfunc(Fai(i,:),fai2.'))^2 + abs(Wfunc(Fai(i,:),fai2'))^2;
        Q(i) = lamda(gg)*U(i,1) + (1-lamda(gg))*U(i,2);
        Q(i) = param.alpha * Q(i) / -2;
    end
    mdl = fitlm(P, Q);
    R2R(gg) = mdl.Rsquared.Adjusted;
    gg
end

% plot(lamda, R2R)
scatter(P(2:num+1), Q(2:num+1),8,[0.4940 0.1840 0.5560], 'filled');  hold on
plot(P,mdl.Coefficients{2,1}*P+mdl.Coefficients{1,1}, 'color', [0.9290 0.6940 0.1250],'linewidth', 1.4);
scatter(P(1), Q(1),8,[0 0.4470 0.7410], 'filled');
scatter(P(2), Q(2),8,[0.8500 0.3250 0.0980], 'filled'); hold off
ylabel('-$\lambda_{R}$H($e^{i\theta}$, $e^{i\varphi_{R1}}$)-$\lambda_{D}$H($e^{i\theta}$, $e^{i\varphi_{D1}}$)', 'FontSize',14,'Interpreter', 'latex');
xlabel('P($\theta$, A)', 'FontSize',14,'Interpreter', 'latex');
title('memorize two patterns', 'FontSize',14,'Interpreter', 'latex');
text(-30, -18, '$r^2=0.947 ^{***}$', 'FontSize',12,'Interpreter', 'latex');

%% 能量函数关系
U = [];
lamda = 0.01:0.01:0.99;
for gg = 1:length(lamda)
    param.Mat = [lamda(gg)*V1 (1-lamda(gg))*V2] * pinv([V1 V2]);  % 连接矩阵
    U(gg,1) = Kuramoto_potential(angle(fai)',param);
    U(gg,2) = Kuramoto_potential(angle(fai2)',param);
end

plot(lamda, U, 'linewidth', 1.4);
xlabel('$\lambda_{R}$ for Weight design', 'Interpreter', 'latex','FontSize',12);
ylabel('P($\theta$, A)', 'Interpreter', 'latex','FontSize',12);
legend('$\theta=\varphi_{R1}$', '$\theta=\varphi_{D1}$','Interpreter', 'latex','FontSize',12);


%% toy landscape
clear;clc
V1 = 0.65;
V2 = 1 - V1;

C1 = [-1 1];
C2 = [1 -1];
C3 = [1 1];
C4 = [-1 -1];

P1 = @(x,y) -V1 * exp( ( (x-C1(1)).^2 + (y-C1(2)).^2 ) ./ -1.4)...
    -V1 * exp( ( (x-C2(1)).^2 + (y-C2(2)).^2 ) ./ -1.4);

P2 = @(x,y) -V2 * exp( ( (x-C3(1)).^2 + (y-C3(2)).^2 ) ./ -1.4)...
    -V2 * exp( ( (x-C4(1)).^2 + (y-C4(2)).^2 ) ./ -1.4);

[X,Y] = meshgrid(-1.3:0.01:1.3, -1.3:0.01:1.3);

A1 = P1(X,Y);

A2 = P2(X,Y);

mesh(A1+A2)

axis off
set(gcf, 'color','none');
