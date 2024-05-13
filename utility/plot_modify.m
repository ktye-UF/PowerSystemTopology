% plot_surrogate_v2(mean_bm, {mean_lhs,mean_pce_exp, mean_pce}, {'LHS','PCE_{ind}','PCE_{cop}'}, 'Mean')
% plot_surrogate_v2(var_bm, {var_lhs,var_pce_exp, var_pce}, {'LHS','PCE_{ind}','PCE_{cop}'}, 'Variance')

% plot_surrogate_v2(var_bm, {var_lhs, var_pce var_pce_exp(1:501)}, {'LHS','PCE_{ind}','PCE_{cop}'}, 'Variance')

% plot_surrogate_v2(var_bm, {var_lhs, var_pce(1,:), var_pce(2,:)}, {'MCS', 'LHS', 'PCE', 'CoPCE'}, 'Variance')

% plot_surrogate_v2(var_bm_out, {var_pce_lhs; var_pce_ind; var_pce_cop;}, {'LHS', 'PCE', 'CoPCE'}, 'Variance')

% close all
ind_machine = 60;
figure; plot(t,x(nge+ind_machine,:)*(180/pi)-x(nge+1,:)*(180/pi),'color','k','linewidth',2.5);
xlabel('Time/s','fontsize',font); 
ylabel(['\delta_{', num2str(ind_machine), '-1} in degree'],'fontsize',font);
title(['Rotor angle \delta_{', num2str(ind_machine), '-1}'])
%% Grid
% % grid
% % gx = geoaxes;
% grid on; grid minor;  % grid minor is a switch?
% % grid off

%% axis
% % label
legend('MCS', 'LHS', 'PCE', 'CoPCE')
% legend('MCS', 'LHS', 'PCE', 'CoPCE')
% legend('MCS', 'LHS', 'MCS', 'MCS')
xlabel('Times (s)', 'FontSize',18, 'FontWeight','bold')
ylabel('Variance of \delta_{40-1} (deg^2)', 'FontSize',18, 'FontWeight','bold')
% ylabel('Mean of \delta_{21} (deg)', 'FontSize',18, 'FontWeight','bold')
% ylabel('Rotor angle \delta_{30-1} (deg)', 'FontSize',18, 'FontWeight','bold')

% % axis line width
set(gca,'linewidth',1.5);  
% % stick font
set(gca,'FontSize',15);

% % axis stick value
xlim([0,10]); ylim([2.5,12]);    % s1_var
xlim([0,10]); ylim([2.5,13.5]); 
xlim([0,10]); ylim([3,12.5]);
xlim([0,10]); ylim([3,14.5]);   % s4
xlim([0,10]); ylim([6,32.5]);     % s5
xlim([0,10]); ylim([0,170]);     % s6: unstable nor
xlim([0,3]); ylim([0,2e4]);     % s6: unstable out
xlim([0,12]); ylim([-200, 0]);     % s8: krig + 6 loads 6 PVs

% % axis minor stick
h = gca; % Minor ticks which don't line up with majors
h.XAxis.MinorTick = 'on'; h.YAxis.MinorTick = 'on'; 
% h.XAxis.MinorTick = 'off'; h.YAxis.MinorTick = 'off'; 
% % h.XAxis.MinorTickValues = 0:1:10; h.YAxis.MinorTickValues = 0:1:32; 


%% figure position
% % figure size
% [left, bottom, width, height]
% set(gcf,'windowstyle','normal');  % for Error: cannot set position when docked
set(gcf,'unit','normalized','position', [0.15,0.15,0.38,0.48]);
% % cover entire window
set(gca,'LooseInset',get(gca,'TightInset'))

error
%% subfigure (manual)
% % off legend and lables
legend off; % axis on % off too mach
% grid off; 
% set(gca,'linewidth',1);   % change all
xlabel(''); ylabel('');     % still there?
set(gca,'xticklabel',''); set(gca,'yticklabel','');

% % off everything
axis off;   % axis on
legend off;


% % xlim: just do it manually
% xlim([4,6]); ylim([2,8])
% xlim([5.1,6.4]); ylim([4,6.5])

%% line color
% lineStyles=linspecer(N,varargin)
lineStyles = linspecer(4)

%% during plot
% % line on top
% h = plot(...)
% uistack(h,'top')  % [error]: will disorder legend



