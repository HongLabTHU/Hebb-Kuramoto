%% memory single phase pattern
clear;clc
node = 36;   % 节点数

F = Traveling_Mode;
grids = F.set_grids(6);

tem1 = F.rotational(grids,[0 0], 1, 1);  % 旋转型
tem2 = F.translational(grids,pi/2,2);  % 前后型

tem1 = conj(tem1');

param.node = node;
endp = 1000;
watchp = 900;


Mat = [real(tem1) imag(tem1)] * pinv([real(tem1) imag(tem1)]);


param.Mat = Mat;  % Connectivity Martix
param.inputtime = 10;    % input time 
param.inputpulse = 0.0;     % input 
param.inputstd = 0.3;
param.sigma = 0.0;
Noise = normrnd(0,2.5,node,endp+2);
Noise(abs(Noise)<1.7) = 0;
param.Noise = [];  % Noise; % OU_noise(1:endp+2, node, 0.1, 200)';
param.a = 1;
param.b = 10 * 2*pi;
param.alpha = 200;
param.dt = 0.01;

init_value = [ones(1,node) 2*pi*rand(1,node)-pi]';
time = 0:param.dt:endp*param.dt;

[t, x] = ode15s(@(t,x) Hopf_fun_tr(t,x, param),time, init_value);
X = x(:, 1 : node).*cos(x(:, node+1 : 2*node));

param.Mat = tem1 * pinv(tem1);
[t, x2] = ode15s(@(t,x) Hopf_fun_tr(t,x, param),time, init_value);

% cm_hsv = colormap('hsv');
% cm_jet = colormap('jet');


%% Phase map 
figure('Position',[488,243.4,320,300])

imagesc(grids(1:6,2), grids(1:6,2),angle(flipud(reshape(tem1,6,6))), [-pi,pi])
hold on
quiver(grids(:,1),grids(:,2),-imag(tem1),-real(tem1), 0.3, 'linewidth',1.5,'color','k');
% quiver(grids(:,1),grids(:,2),zeros(node,1),-ones(node,1), 0.3, 'linewidth',1.5,'color','k');
colorbar('Ticks',[-pi, 0, pi], 'TickLabels',{'-\pi','0','\pi'}, 'FontSize',12);
grid on;axis off
title('R1', 'FontSize',12, 'Interpreter', 'latex');
colormap hsv

%% Phase map 
figure
imagesc(1:200,1:node,angle(exp(1j*x(1:200,node+1:2*node)')),[-pi,pi])
colormap(gca,winter);
colorbar('Ticks',[-pi, 0, pi], 'TickLabels',{'-\pi','0','\pi'}, 'FontSize',12);
xlabel('timepoints', 'FontSize',12, 'Interpreter', 'latex');
ylabel('nodes', 'FontSize',12, 'Interpreter', 'latex');
title('Connectivity strength based method', 'FontSize',12, 'Interpreter', 'latex');

fai3 = x2(watchp, node+1:2*node);
fai3 = exp(1i*fai3);
% 'delay'
abs(Wfunc(fai3,conj(tem1')))
fai3 = x(watchp, node+1:2*node);
fai3 = exp(1i*fai3);
% 'Connectivity'
abs(Wfunc(fai3,conj(tem1')))


%% Phase 
figure
imagesc(0.5:1:node+0.5, 0.5:1:node+0.5, Mat);
% colormap(gca,hot);
colorbar
title('Connectivity matrix', 'FontSize',12, 'Interpreter', 'latex');
grid on

figure
polarscatter(angle(fai3), abs(fai3), 14, [0 0 0], 'filled')
set(gca,'rTicklabel',[])
thetaticks(gca,[0 45 90 135 180 225 270 315])
thetaticklabels(gca,{'0','\pi/4','\pi/2','3\pi/4','2\pi','-3\pi/4','-\pi/2','-\pi/4'});
title('Memorized pattern', 'FontSize',12, 'Interpreter', 'latex');

polarscatter(angle(Mat*conj(fai3')), abs(Mat*conj(fai3')), 14, [0 0 0], 'filled')
title('Endpoint pattern', 'FontSize',12, 'Interpreter', 'latex');


cmap = colormap('hsv');
for i = 0.5:0.01:1
    for j = -pi:0.1:pi
        scatter(i * cos(j), i* sin(j), [], cmap(floor((j+pi)/(2*pi) * 255 + 1), :), 'filled');
        hold on
    end
end

%% Potential function
dt = param.dt;
P = zeros(2, endp);
for i = 1:endp
%     P(1, i) = Kuramoto_potential(x(i ,node+1:2*node), param);
%     P(2, i) = Kuramoto_potential(x2(i ,node+1:2*node), param);
    P(1, i) = abs(Wfunc(exp(1j*x(i ,node+1:2*node)),conj(tem1')));
    P(2, i) = abs(Wfunc(exp(1j*x2(i ,node+1:2*node)),conj(tem1')));
end

% P1 = Kuramoto_potential(angle(tem1), param)*ones(1,length(P));

plot(P','linewidth',1.2);
xlabel('timepoints', 'FontSize',12, 'Interpreter', 'latex');
legend('Connectivity strength based method','Connectivity delay based method',  'FontSize',10,'Interpreter', 'latex');
xlabel('timepoints', 'FontSize',12, 'Interpreter', 'latex');
ylabel('$$ \left| \rho(exp(i\theta), \Phi_{R1}) \right| $$', 'FontSize',12, 'Interpreter', 'latex');