%% �켣
clear;clc
node = 36;   % �ڵ���

F = Traveling_Mode;
grids = F.set_grids(6);

tem1 = F.rotational(grids,[0 0], 1, 1);  % ��ת��
tem2 = F.translational(grids,pi/2,2);  % ǰ����

tem = [tem1; conj(tem1); tem2; conj(tem2)];
tem = tem.';

param.node = node;
endp = 1200;
watchp = 900;

fai = tem(:,1);
fai2 = tem(:,3);
V1 = [real(fai) imag(fai)];
V2 = [real(fai2) imag(fai2)];
Mat = [0.6*V1 0.3*V2] * pinv([V1 V2]);  % ���Ӿ���


param.Mat = Mat;  % ���Ӿ���
param.inputtime = 10;    % ����ʱ��
param.inputpulse = 0.0;     % �����С
param.inputstd = 0.3;
param.sigma = 0.0;
Noise = normrnd(0,2.5,node,endp+2);
Noise(abs(Noise)<1.7) = 0;
param.Noise = [];  % Noise; % OU_noise(1:endp+2, node, 0.1, 200)';
param.a = 1;
param.b = 10 * 2*pi;
param.alpha = 200;
param.dt = 0.01;

time = 0:param.dt:endp*param.dt;
X = zeros(100, endp+1, node);

for bash = 1:100
    init_value = [ones(1,node) 2*pi*rand(1,node)-pi]';
    [t, x] = ode15s(@(t,x) Hopf_fun(t,x, param),time, init_value);
    X(bash,:,:) = x(:, node+1 : 2*node);
end

P = zeros(100, endp);
for bash = 1:100
    for j = 2:(endp+1)
        P(bash, j-1) = Kuramoto_potential(squeeze(X(bash, j ,:)), param);
    end
end

figure('Position',[766.6,370.6,296,269.6])
for bash = 1:100
    v = squeeze(exp( 1i * X(bash,2:end,:) ));
    res = wocc_projection(v, tem, 2);
    % Colr = [160*360/239, 0/240, 197/240];
    % Colr = [28*360/239, 204/240, 127/240];
%     if norm(res(end, :) - [1 0]) < 2e-1
%         Colr = [136*360/239, 240/240, 89/240];
%     end
%     if norm(res(end, :) - [-1 0]) < 2e-1
%         Colr = [12*360/239, 190/240, 114/240];
%     end
%     if norm(res(end, :) - [0 1]) < 2e-1
%         Colr = [28*360/239, 204/240, 127/240];
%     end
%     if norm(res(end, :) - [0 -1]) < 2e-1
%         Colr = [57*360/239, 135/240, 104/240];
%     end
%     Colr = [Colr(1)*ones(endp,1) Colr(2)*ones(endp,1) ...
%         linspace(0.8+0.2*Colr(3), 0.3+0.7*Colr(3), endp)'];
    
%     index = floor(endp*(max(P(bash,:))-P(bash,:)) / (max(P(bash,:)) - min(P(bash,:))))+1;
%     index(index>endp) = endp;
%     [R_, G_, B_] = hsl2rgb(Colr(index,1), Colr(index,2), Colr(index,3));
%     scatter(res(:,1), res(:,2), 0.6, [R_, G_, B_]);
    scatter(res(:,1), res(:,2), 0.6, log10(1:endp));
    colormap summer
    hold on
end


axis([-1 1 -1 1])
xlabel('W(exp(i$\theta)$, $ \Phi_{R})$', 'FontSize',12, 'Interpreter', 'latex');
ylabel('W(exp(i$\theta)$, $ \Phi_{D})$', 'FontSize',12, 'Interpreter', 'latex');
title('$\lambda_{R}$=0.6, $\lambda_{D}$=0.3', 'FontSize',12, 'Interpreter', 'latex');

%% �ֲ�ͼ
Lamd = 0.01:0.001:0.99;
G = zeros(72,length(Lamd));
H = zeros(1000,length(Lamd));
i = 1;
for lamd = Lamd
    for alp = 200
        Mat = [lamd*V1 (1-lamd)*V2] * pinv([V1 V2]);  % ���Ӿ���
        param.Mat = Mat;
        param.alpha = alp;
        param.node = 36;
        J1 = Kuramoto_Jakobi(angle(fai),param, 2);
        J1 = real(eig(J1));
        J1 = sort(J1,'descend');
        % J1(J1==0)=[];
        J2 = Kuramoto_Jakobi(angle(fai2),param, 2);
        J2 = real(eig(J2));
        J2 = sort(J2,'descend');
        % J2(J2==0)=0;
        G(:,i) = [J1; J2];
        H(alp, i) = J1(1);
        H(alp + 500, i) = J2(1);
    end
    i = i +1 
end

colormap jet;
imagesc(Lamd, 0.5:1:35.5, G(1:36,:), [-6 6]);
set(gca,'ygrid','on')
set(gca, 'YTick', 1:36, 'YTickLabel', []);
xlabel('$\lambda_{R}$ for Weight design',  'FontSize',12, 'Interpreter', 'latex');
ylabel('eigenvalue of $\Phi_{R}$',  'FontSize',12, 'Interpreter', 'latex');
xline(Lamd(385), '--m', 'linewidth', 1);
% xline(Lamd(440), '--', 'color', [0.9290 0.6940 0.1250], 'linewidth', 1);


imagesc(Lamd, 0.5:1:35.5, G(37:72,:), [-6 6]);
set(gca,'ygrid','on')
set(gca, 'YTick', 1:36, 'YTickLabel', []);
xlabel('$\lambda_{R}$ for Weight design', 'FontSize',12, 'Interpreter', 'latex');
ylabel('eigenvalue of $\Phi_{D}$', 'FontSize',12,  'Interpreter', 'latex');
% xline(Lamd(295), '--m', 'linewidth', 1);
xline(Lamd(530), '--', 'color', [0.9290 0.6940 0.1250], 'linewidth', 1);
imagesc(G,[-6 6])
colorbar('Ticks',[-6, -3, 0, 3, 6], 'TickLabels',{'-6','-3', '0','3', '6'});


%% Բ�ܲ�����
Phi = 0:pi/200:2*pi;
X = imag(tem(:,3)); Y = real(tem(:,1)); Z = real(tem(:,3));
% scatter(grids(:,1),grids(:,2),[],Z);
InitD = zeros(length(Phi), node);
i = 1;
for phi = Phi
    % פ��
    U = X;
    V = (cos(phi)*Y + sin(phi)*Z);
    InitD(i, :) = exp(1j*angle(V+1j*U)).';
    %     V = (-sin(phi)*Y + cos(phi)*Z);
    %     InitD(i, :) = exp(1j*angle(V+1j*U)).';
    res = wocc_projection(InitD(i, :), tem, 2);
    K = Kuramoto_potential(angle(InitD(i, :)), param);
    i= i + 1;
    hold on
    scatter(res(1), res(2), 15, K, 'filled');
end

xlabel('W(exp(i$\theta)$, exp(i$\varphi_{R})$)', 'FontSize',12, 'Interpreter', 'latex');
ylabel('W(exp(i$\theta)$, exp(i$\varphi_{D})$)', 'FontSize',12, 'Interpreter', 'latex');
title('$\lambda_{R}$=0.59, $\lambda_{D}$=0.41', 'FontSize',12, 'Interpreter', 'latex');

%%
K = zeros(3, size(InitD, 1));
for i = 1:size(InitD, 1)
    param.Mat = [0.36*V1 0.64*V2] * pinv([V1 V2]);  % ���Ӿ���
    K(1, i) = Kuramoto_potential(angle(InitD(i, :)), param);
    param.Mat = [0.48*V1 0.52*V2] * pinv([V1 V2]);  % ���Ӿ���
    K(2, i) = Kuramoto_potential(angle(InitD(i, :)), param);
    param.Mat = [0.59*V1 0.41*V2] * pinv([V1 V2]);  % ���Ӿ���
    K(3, i) = Kuramoto_potential(angle(InitD(i, :)), param);
end
plot(K')
