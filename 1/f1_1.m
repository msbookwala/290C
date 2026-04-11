% MAE 190: CFD, Programming Assignment #1
% OTS, 2020

clear variables
close all
clc

% load data from file
load('cylinder_Re100.mat')

% determine number of snapshots and grid size
[nt,nx,ny]  = size(u);
% NOTE: the x- and y- indices correspond to Matlab's 'ndgrid' format.
% ndgrids can be much nicer to work with, as it sticks with conventional
% indexing. That is, the first index corresponds to x, and the second to y,
% etc. This is different from the'meshgrid' format. It can get a bit confusing
% between the two, so just be sure to keep track of what kind of grid
% is being used.

%% 1.1: Basic flow visualization 
figure
for ti = 1:nt
    % plot u
    subplot(2,1,1)
    pcolor(x,y,squeeze(u(ti,:,:)))
    shading interp, axis equal tight
    caxis([-0.5 2])
    title(['u (' num2str(ti) '/' num2str(nt) ')']); xlabel('x'); ylabel('y')
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off
    %colorbar
    
    % plot v
    subplot(2,1,2)
    pcolor(x,y,squeeze(v(ti,:,:)))
    shading interp, axis equal tight 
    caxis([-1 1])
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off
    title(['v (' num2str(ti) '/' num2str(nt) ')']); xlabel('x'); ylabel('y')
    %colorbar
    drawnow
end

%% 1.2 mean
% Calculate mean velocity fields
meanU = squeeze(mean(u(151:end,:,:), 1));
meanV = squeeze(mean(v(151:end,:,:), 1));

% Plot mean velocity fields
figure
subplot(2,1,1)
pcolor(x, y, squeeze(meanU));
shading interp, axis equal tight
title('Mean U Velocity'); xlabel('x'); ylabel('y')
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off
subplot(2,1,2)
pcolor(x, y, squeeze(meanV));
shading interp, axis equal tight
title('Mean V Velocity'); xlabel('x'); ylabel('y')
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off
%% 1.3 streamline

% Plot mean velocity fields
figure
pcolor(x, y, squeeze(meanU));
shading interp, axis equal tight
title('Mean U Velocity'); xlabel('x'); ylabel('y')
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1])
streamline(x',y',meanU',meanV',-5*ones(10,1),[-4 -3 -2 -1 -0.01 0.01 1 2 3 4]); hold off

%% 1.4 

figure
for ti = nt:nt
    % plot u
    subplot(2,1,1)
    pcolor(x,y,squeeze(u(ti,:,:))-meanU)
    shading interp, axis equal tight
    caxis([-0.5 0.5])
    title(['u (' num2str(ti) '/' num2str(nt) ')']); xlabel('x'); ylabel('y')
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off
    %colorbar
    
    % plot v
    subplot(2,1,2)
    pcolor(x,y,squeeze(v(ti,:,:))-meanV)
    shading interp, axis equal tight 
    caxis([-0.5 0.5])
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off
    title(['v (' num2str(ti) '/' num2str(nt) ')']); xlabel('x'); ylabel('y')
    %colorbar
    drawnow
end
%% 1.5
meanU_ = permute(meanU,[3,1,2]);
meanV_ = permute(meanV,[3,1,2]);

k = 0.5*(mean((u(151:end,:,:)-meanU_).^2, 1)+ mean((v(151:end,:,:)-meanV_).^2, 1));
pcolor(x,y,squeeze(k))
shading interp, axis equal tight 
rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]); hold off
title('Turbulent Kinetic Energy');xlabel('x');ylabel('y')
drawnow