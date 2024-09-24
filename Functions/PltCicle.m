function [x,y] = PltCicle(center,radius)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% 定义圆心和半径
% center = [2, 3]; % 圆心坐标
% radius = 5;      % 圆半径

% 生成角度数组
theta = linspace(0, 2*pi, 100);

% 计算圆上的点的坐标
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);


% % 画圆
% figure;
% plot(x, y, 'r', 'LineWidth', 2);
% hold on;
% 
% % 标记圆心
% plot(center(1), center(2), 'ko', 'MarkerSize', 8);

end