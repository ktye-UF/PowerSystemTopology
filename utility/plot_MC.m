function [] = plot_MC(Y_bm, num_sample)
global t_end ind_output
mean_bm = mean(Y_bm); var_bm = var(Y_bm);
% % 
figure; hold on
t = 0:0.02:t_end;

plot(t, mean_bm-3.*sqrt(var_bm), '--b', 'linewidth', 2);
plot(t, mean_bm+3.*sqrt(var_bm), '--r', 'linewidth', 2);

% plot([t(1), t(end)], [180, 180], '--');
% plot([t(1), t(end)], [-180, -180], '--')
%%
% % number of sample
% for p=1:num_sample
%     plot(t, Y_bm(p,:), 'g'); 
% %     plot(t, Y_bm(p,:)); 
% end
% % rand
if num_sample > size(Y_bm,1)
    plot(t, Y_bm, 'g')
else
    plot(t, Y_bm(randperm(size(Y_bm,1),num_sample),:), 'g'); 
end


%%

if length(ind_output) == 2
    ylabel(['Rotor angle \delta_{', num2str(ind_output(1)), '-', num2str(ind_output(2)), '} (deg)'])
%     title(['mean of \delta_{', num2str(ind_output(1)), num2str(ind_output(2)), '}'])
else
    ylabel(['Rotor speed \omega_{', num2str(ind_output), '} (pu)'])
%     title(['mean of \omega_{', num2str(ind_output), '}'])
end
xlabel('Times (s)', 'FontSize',18, 'FontWeight','bold')
ylabel('Rotor angle of \delta_{40-1} (deg)', 'FontSize',18, 'FontWeight','bold')

set(gca,'linewidth',1.5); 
set(gca,'FontSize',15);
xlim([0,6]); %ylim([2.5,12]);    % s1_var
h = gca; % Minor ticks which don't line up with majors
h.XAxis.MinorTick = 'on'; h.YAxis.MinorTick = 'on'; 
set(gcf,'unit','normalized','position', [0.15,0.15,0.38,0.48]);
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca, 'child', [h3,h1,h2])
legend('Lower bound', 'Upper bound', 'MCS')

plot(t, mean_bm, '--k', 'linewidth', 1.5); 




% % 
% legend('Upper bound', 'Lower bound', 'Mean', 'MCS')
% legend('Mean', 'Upper bound', 'Lower bound', 'MCS')
% uistack(h,'top')
% ch=get(gca,'children');
% set(gca,'children',[ch(2) ch(1)]);

% ylabel('degree')
% ylim([-200, 200])
% title('correlated PV [0.8,0.6,0.4]'); title('correlated PV [0,0,0]')
% 
% figure; errorbar(t,mean_bm,sqrt(var_bm))
% figure; plot(t, var_bm)










