%% macroscopic recordings (LFP)---Figure 2

figure('Position',[0 0 1200 800]);
% Figure 2a
subplot(311)
plot(t+t0,multiplyer*Ve_mea_nd(1,:),Color='b',LineWidth=1); 
hold on
plot(t+t0,multiplyer*Ve_mea_nd(2,:),'Color',[0, 0.8, 0],LineWidth=1); 
hold on
plot(t+t0,multiplyer*Ve_mea_nd(3,:),Color='r',LineWidth=1); 
hold on
xlabel('time (ms)')
ylabel('LFP (mV)')
ylim([-0.1 0.6])
leg=legend('V_A','V_B','V_C');
leg.FontWeight = 'bold';
leg.FontSize = 12;
gca_width=1.2;
set(gca,'linewidth',gca_width)
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize=15;
ax.YAxis.FontSize=15;
title('Figure 2a')
% Figure 2b
subplot(312)
gca_width=1.2;
plot(t,Id2_in,LineWidth=1.2)
ylim([0 20])
xlabel('time (ms)')
ylabel('I_{d2} (mV)')
set(gca,'linewidth',gca_width)
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize=15;
ax.YAxis.FontSize=15;
title('Figure 2b')
% Figure 2d
Fs = 20000;
signal_1=Ve_mea_nd(1,:);
signal_2=Ve_mea_nd(2,:);
signal_3=Ve_mea_nd(3,:);
fc = 1;                 
order = 4;              
Wn = fc/(Fs/2); 
[b, a] = butter(order, Wn, 'low');
filtered_signal_1 = filtfilt(b, a, signal_1);
filtered_signal_2 = filtfilt(b, a, signal_2);
filtered_signal_3 = filtfilt(b, a, signal_3);
subplot(313)
plot(t+t0,filtered_signal_1,Color='b',LineWidth=1.5); 
hold on
plot(t+t0,filtered_signal_2,'Color',[0, 0.8, 0],LineWidth=1.5); 
hold on
plot(t+t0,filtered_signal_3,Color='r',LineWidth=1.5); 
hold on
xlabel('time (ms)')
ylabel('LFP (mV)')
ylim([-0.15 0.15])
leg=legend('V_A','V_B','V_C');
leg.FontWeight = 'bold';
leg.FontSize = 12;
gca_width=1.2;
set(gca,'linewidth',gca_width)
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize=15;
ax.YAxis.FontSize=15;
title('Figure 2d')




%% microscopic recordings (membrane potential and endogenous electric field)---Figure 3

% Figure 3a
gca_width=1.2;
x=[1:nN];
ts=2000;te=2500;
[X,Y]=meshgrid(t(ts/dt:te/dt),x);
figure('Position', [0, 0, 400, 600]);
subplot(211)
surf(X,Y,Vmd1(:,ts/dt:te/dt));
shading flat
view(2);
cb=colorbar;
set(cb, 'FontSize', 10, 'FontWeight', 'bold');
colormap(jet);
xlim([ts, te])
xlabel("time (ms)")
ylabel("Cell number")
set(gca,'YDir','normal')
set(gca,'linewidth',gca_width)
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize=10;
ax.YAxis.FontSize=10;
title('Figure 3a1')
subplot(212)
surf(X,Y,Vmd2(:,ts/dt:te/dt));
shading flat
view(2);
cb=colorbar;
set(cb, 'FontSize', 10, 'FontWeight', 'bold');
colormap(jet);
xlim([ts, te])
xlabel("time (ms)")
ylabel("Cell number")
set(gca,'YDir','normal')
set(gca,'linewidth',gca_width)
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize=10;
ax.YAxis.FontSize=10;
title('Figure 3a2')

% Figure 3c
figure('Position', [0, 0, 400, 600]);
subplot(211)
surf(X,Y,E_ext_sd1(:,ts/dt:te/dt));
shading flat
view(2);
cb=colorbar;
set(cb, 'FontSize', 10, 'FontWeight', 'bold');
colormap(jet);
xlim([ts, te])
xlabel("time (ms)")
ylabel("Cell number")
set(gca,'YDir','normal')
set(gca,'linewidth',gca_width)
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize=10;
ax.YAxis.FontSize=10;
title('Figure 3c1')
subplot(212)
surf(X,Y,E_ext_sd2(:,ts/dt:te/dt));
shading flat
view(2);
cb=colorbar;
set(cb, 'FontSize', 10, 'FontWeight', 'bold');
colormap(jet);
xlim([ts, te])
xlabel("time (ms)")
ylabel("Cell number")
set(gca,'YDir','normal')
set(gca,'linewidth',gca_width)
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize=10;
ax.YAxis.FontSize=10;
title('Figure 3c2')












