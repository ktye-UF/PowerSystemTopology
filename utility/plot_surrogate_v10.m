function [] = plot_surrogate_v10(x, y, sur_name)
%% v10 -> for ppf (instead of dynamic)
% % stat_bm: N*T 
% % stat_sur: {N*T, N*T, ...}
% % sur_name: {}
% % stat_name: angle or speed
% global t_end ind_gen
% t = 0:0.02:t_end;
n_curve = size(x,1);
% [MCS, LHS, PCE, CoPCE]
% plot_line = {':', '--', '-.'};
plot_line = {'-', '-', ':', '--', '-.'};
line_width = [3, 2.5, 2.4, 2.2, 2.5];
% line_width = [2.5, 2, 2.2, 1.5, 2];
% % red color
% line_color = [0.88,0.63,0.62]; [0.46,0.14,0.14]; [0.32,0.19,0.19]; [0.63,0.18,0.18];  
% % multi color
line_color = linspecer(n_curve);
line_color([end-1,end],:) = line_color([end,end-1],:);
% line_color = {[0.88,0.63,0.62], [0.46,0.14,0.14], [0.32,0.19,0.19], [0.63,0.18,0.18]};  

%% 
figure; hold on
str_leg = {};
for i = 1:n_curve
    plot(x(i,:), y(i,:), plot_line{i}, 'color', line_color(i,:), 'linewidth', line_width(i))
    str_leg = [str_leg, sur_name{i}];
end
legend(str_leg);
xlabel('Line Power Flow (MVA)', 'FontSize',18, 'FontWeight','bold')
ylabel('Probability Density Function', 'FontSize',18, 'FontWeight','bold')

%%
% figure;hold on
% plot(t, stat_bm, plot_line{1}, 'color', line_color(1,:), 'linewidth', line_width(1))
% str_leg = {'MCS'};
% for i=1:n_sur
% %     plot(t, stat_sur{i}, [plot_line{i},'b'])
%     plot(t, stat_sur(i,:), plot_line{i+1}, 'color', line_color(i+1,:), 'linewidth', line_width(i+1))
% %     plot(stat_sur{i}, plot_line{i+1}, 'color', line_color(i+1,:), 'linewidth', line_width(i+1))
%     str_leg = [str_leg, sur_name{i}];
% end
% legend(str_leg);
% xlabel('Times (s)', 'FontSize',18, 'FontWeight','bold')
% ylabel('Variance of rotor angle \delta_{40-1} (deg^2)', 'FontSize',18, 'FontWeight','bold')

% ylabel([stat_name, ' of \delta_{', num2str(ind_gen), '-1}'])
% if length(ind_output) == 2
%     ylabel([stat_name, ' of \delta_{', num2str(ind_output(1)), num2str(ind_output(2)), '}'])
% else
%     ylabel([stat_name, ' of \omega_{', num2str(ind_output), '} (pu)'])
% end
set(gca,'linewidth',1.5); 
set(gca,'FontSize',15);
% xlim([0,6]); %ylim([2.5,12]);    % s1_var
h = gca; % Minor ticks which don't line up with majors
h.XAxis.MinorTick = 'on'; h.YAxis.MinorTick = 'on'; 
set(gcf,'unit','normalized','position', [0.15,0.15,0.38,0.48]);
set(gca,'LooseInset',get(gca,'TightInset'))







