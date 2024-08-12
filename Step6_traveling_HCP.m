%% HCP 复现CPCA行波
% 查看原始行波
clear;clc;
load('.\utility\CPCAsurf.mat');
for i = 1:3
    figure
    trimesh(Face_L+1., Vert_L(:,2),Vert_L(:,1),Vert_L(:,3), pca_comps(i,1:2562),'LineWidth',1.5);
    colormap hsv
    caxis([-pi,pi])
    title(['CPCA comps ' num2str(i)])
end
colorbar


%% 参数设置
node = size(pca_comps,2);
V1 = [cos(pca_comps(1,:)') sin(pca_comps(1,:)')];
V2 = [cos(pca_comps(2,:)') sin(pca_comps(2,:)')];
V3 = [cos(pca_comps(3,:)') sin(pca_comps(3,:)')];

param.node = node;
endp = 30;
watchp = 25;

Mat = [V1 V2 V3] * pinv([V1 V2 V3]);

param.Mat = Mat;  % 连接矩阵
param.inputtime = 10;    % 输入时间
param.inputpulse = 0.0;     % 输入大小
param.inputstd = 0.3;
param.sigma = 0.0;

param.Noise = [];
param.a = 1;
param.b = 10 * 2*pi;
param.alpha = 8*node;
param.dt = 0.01;

% 欧拉法
Af = zeros(1000, 3);
x = 2*pi*rand(1,node)'-pi;
X = zeros(node,1002);
U = sin(x - x');
tic

for t = 1:1001
    U = sin(x - x');
    dx = param.alpha * diag(Mat * U) / node;
    x = x + 0.002 * dx;
    X(:,t+1) = x;
    Af(t,1) = abs(Wfunc(exp(1j*x'), exp(1j*pca_comps(1,:))))^2 -...
        abs( Wfunc( exp(1j*x'), exp(-1j*pca_comps(1,:)) ) )^2;
    Af(t,2) = abs(Wfunc(exp(1j*x'), exp(1j*pca_comps(2,:))))^2 -...
        abs( Wfunc( exp(1j*x'), exp(-1j*pca_comps(2,:)) ) )^2;
    Af(t,3) = abs(Wfunc(exp(1j*x'), exp(1j*pca_comps(3,:))))^2 -...
        abs( Wfunc( exp(1j*x'), exp(-1j*pca_comps(3,:)) ) )^2;
    t
end

toc


%% 历时 14125.090263 秒。
tic
inT = 0.23:0.02:0.43;
B = zeros(length(inT),length(inT));
for x = 1:length(inT)
    for y = 1:length(inT)
        for i = 1:3
            param.Mat = [inT(x)*V1 inT(y)*V2 (1-inT(x)-inT(y))*V3] * pinv([V1 V2 V3]);
            J = Kuramoto_Jakobi(pca_comps(i,:),param, 2);
            D = eig(J);
            if max(real(D)) > 2e-10
                B(x,y) = B(x,y) + 2^(i-1);
            end
%             figure
%             scatter(real(D),imag(D),'filled');
%             xline(0);
%             title(num2str(i))
        end
        toc
    end
    x
end

imagesc(inT,inT,flipud(B))
set(gca,'YTick', 0.24:0.02:0.42);
set(gca,'YTicklabel', {'0.42','0.40','0.38','0.36','0.34',...
    '0.32','0.30', '0.28', '0.26', '0.24'});
% set(gca,'XTicklabel', {'0.22','0.24','0.26','0.28','0.30','0.32',...
%     '0.34','0.36', '0.38', '0.40', '0.42', '0.44'});
% xlabel('$\lambda_{P2}$', 'FontSize',12, 'Interpreter', 'latex');
% ylabel('$\lambda_{P1}$', 'FontSize',12, 'Interpreter', 'latex');
% grid on
% 000 0
% 001 1
% 010 2
% 011 3
% 100 4
% 101 5
% 110 6

%% GPU 加速版
% pca_comps = gpuArray(pca_comps);
% pca_comps = single(pca_comps);
node = size(pca_comps,2);
V1 = [cos(pca_comps(1,:)') sin(pca_comps(1,:)')];
V2 = [cos(pca_comps(2,:)') sin(pca_comps(2,:)')];
V3 = [cos(pca_comps(3,:)') sin(pca_comps(3,:)')];

tic
inT = 0.23:0.01:0.43;
B = zeros(length(inT),length(inT));
for x = 1:length(inT)
    for y = 1:length(inT)
        for i = 1:3
            M = [inT(x)*V1 inT(y)*V2 (1-inT(x)-inT(y))*V3] * pinv([V1 V2 V3]);
            Q = M .* cos(pca_comps(i,:) - pca_comps(i,:)');
            J = Q - diag(sum(Q,2));
            D = eig(J);
            if max(real(D)) > 2e-10
                B(x,y) = B(x,y) + 2^(i-1);
            end
%             figure
%             scatter(real(D),imag(D),'filled');
%             xline(0);
%             title(num2str(i))
        end
        toc
    end
    x
end

colormap pink
s = pcolor(B);
s.EdgeColor = [0.8 0.8 0.8];
set(gca,'XTick', 1.5:2:21.5);
set(gca,'XTicklabel', {'0.23','0.25','0.27','0.29','0.31','0.33',...
   '0.35','0.37', '0.39', '0.41', '0.43'});
set(gca,'YTick', 1.5:2:21.5);
set(gca,'YTicklabel', {'0.23','0.25','0.27','0.29','0.31','0.33',...
   '0.35','0.37', '0.39', '0.41', '0.43'});
xlabel('$\lambda_{P2}$', 'FontSize',16, 'Interpreter', 'latex');
ylabel('$\lambda_{P1}$', 'FontSize',16, 'Interpreter', 'latex');

%%
hold on
B_ = B';
for i = 2:20
    for j = 1:20
    if (B_(i,j) - B_(i-1, j) ~=0 )
        if (ismember(B_(i,j), [0 2 4 6]) && ~ismember(B_(i-1, j), [0 2 4 6])) || ...
            (~ismember(B_(i,j), [0 2 4 6]) && ismember(B_(i-1, j), [0 2 4 6]))
            line([i i], [j j+1], 'color', [0 0.4470 0.7410], 'LineWidth', 1.5);
        end
        if (ismember(B_(i,j), [0 1 4 5]) && ~ismember(B_(i-1, j), [0 1 4 5])) || ...
            (~ismember(B_(i,j), [0 1 4 5]) && ismember(B_(i-1, j), [0 1 4 5]))
            line([i i], [j j+1], 'color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
        end
        if (ismember(B_(i,j), [0 1 2 3]) && ~ismember(B_(i-1, j), [0 1 2 3])) || ...
            (~ismember(B_(i,j), [0 1 2 3]) && ismember(B_(i-1, j), [0 1 2 3]))
            line([i i], [j j+1], 'color',[0.9290 0.6940 0.1250], 'LineWidth', 1.5);
        end
    end
    end
end


for i = 1:20
    for j = 2:20
    if (B_(i,j) - B_(i, j-1) ~=0 )
        if (ismember(B_(i,j), [0 2 4 6]) && ~ismember(B_(i, j-1), [0 2 4 6])) || ...
            (~ismember(B_(i,j), [0 2 4 6]) && ismember(B_(i, j-1), [0 2 4 6]))
            line([i i+1], [j j], 'color',[0 0.4470 0.7410], 'LineWidth', 1.5);
        end
        if (ismember(B_(i,j), [0 1 4 5]) && ~ismember(B_(i, j-1), [0 1 4 5])) || ...
            (~ismember(B_(i,j), [0 1 4 5]) && ismember(B_(i, j-1), [0 1 4 5]))
            line([i i+1], [j j], 'color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
        end
        if (ismember(B_(i,j), [0 1 2 3]) && ~ismember(B_(i, j-1), [0 1 2 3])) || ...
            (~ismember(B_(i,j), [0 1 2 3]) && ismember(B_(i, j-1), [0 1 2 3]))
            line([i i+1], [j j], 'color',[0.9290 0.6940 0.1250], 'LineWidth', 1.5);
        end
    end
    end
end




figure
h=trimesh(Face_L+1., Vert_L(:,3),Vert_L(:,1),Vert_L(:,2), angle(y(1415,1:2562)),'LineWidth',1.5);
colormap hsv
caxis([-pi,pi])
set(h, 'FaceAlpha', 0.1);  % 将表面的透明度设置为0.5
axis equal;  % 设置坐标轴比例一致，使球体看起来更圆滑
grid off;
view(0, 90);
