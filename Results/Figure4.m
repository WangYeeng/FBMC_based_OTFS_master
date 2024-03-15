%==========================================================================
%Fig.4
%==========================================================================
clear; close all;clc;
%% Import data
load ('fNormalized.mat');
load ('PSD_FBMC.mat');
load ('PSD_OFDM.mat');
load ('PSD_OFDM_OTFS.mat');
load ('PSD_FBMC_OTFS_Kaiser0.mat');
PSD_FBMC_OTFS_Kaiser0 = PSD_FBMC_OTFS;
load ('PSD_FBMC_OTFS_Kaiser4.mat');
PSD_FBMC_OTFS_Kaiser4 = PSD_FBMC_OTFS;
load ('PSD_FBMC_OTFS_Kaiser9.mat');
PSD_FBMC_OTFS_Kaiser9 = PSD_FBMC_OTFS;
load ('PSD_FBMC_OTFS_Nuttall.mat');
PSD_FBMC_OTFS_Nuttall = PSD_FBMC_OTFS;
load ('PSD_FBMC_OTFS_Tukey1.mat');
PSD_FBMC_OTFS_Tukey1 = PSD_FBMC_OTFS;

%% Plot Results
LineWidth = 1.2;
figure(2)
plot(f_Normalized,PSD_OFDM_OTFS ,'Color',0.80*[1,0,0],'LineWidth',LineWidth);hold on;grid on;
% plot(f_Normalized,PSD_OFDM ,     'Color',0.85*[0,0,0],'LineWidth',LineWidth);
plot(f_Normalized,PSD_FBMC ,     'Color',0.75*[0,1,0],'LineWidth',LineWidth);
plot(f_Normalized,PSD_FBMC_OTFS_Kaiser0 ,'Color',0.75*[1,0,1], 'LineWidth',LineWidth);
plot(f_Normalized,PSD_FBMC_OTFS_Kaiser4 ,'Color',0.75*[1,1,0], 'LineWidth',LineWidth);
plot(f_Normalized,PSD_FBMC_OTFS_Kaiser9 ,'Color',0.75*[0,0,1], 'LineWidth',LineWidth);
%plot(f_Normalized,PSD_FBMC_OTFS_Nuttall ,'Color',0.85*[0,0,1], 'LineWidth',LineWidth);
%plot(f_Normalized,PSD_FBMC_OTFS_Tukey1  ,'Color',0.85*[0,0,1], 'LineWidth',LineWidth);
ylim([-150 5]);
xlim([-5 15]);
xlabel('Normalized frequency $f/F$ (Hz)','Interpreter','latex');
ylabel('Power spectral density (dB)');     
set(gca,'YTick',[-200 -150 -100 -50 -0],'FontName','Times New Roman','FontSize',12,'GridLineStyle','-.');

