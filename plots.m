close all
clear
load('build\bin\values.mat')

SUBSET = 5001;
ENV = 0.6;
FONTSIZE = 12;
LINEWIDTH = 1.5;

% Changing the size of the legend line is overwritten by the line width of the plot
% Don't know if it's MATLAB online or it's 20b (it was working on 19b)

%%
figure
hold on

plot(q_dot_ref(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(qs_dot(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(qm_dot(1:SUBSET), 'LineWidth', LINEWIDTH);
[~, icons] = legend('$\dot{q}_{ref}$','$\dot{q}_{s}$','$\dot{q}_{m}$','Interpreter','latex','FontSize', FONTSIZE);
for idx=1:9
    icons(idx).LineWidth = 6;
end

grid on
title('Comparison between reference and actual trajectories')
hold off

% print('-dpng', '-r300', 'vel')
%%
figure
hold on

plot(q_ref(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(qs(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(qm(1:SUBSET), 'LineWidth', LINEWIDTH);
plot([0,SUBSET],[ENV,ENV], 'LineWidth', LINEWIDTH);
[~, icons] = legend('$q_{ref}$','$q_{s}$','$q_{m}$','env','Interpreter','latex','FontSize', FONTSIZE);
for idx=1:12
    icons(idx).LineWidth = 6;
end

grid on
title('Comparison between reference and actual trajectories')
hold off

% print('-dpng', '-r300', 'pos')
%%
figure 
hold on

plot(tau_plm(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(tau_tlm(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(tau_tlc(1:SUBSET), 'LineWidth', LINEWIDTH);
[~, icons] = legend('$\tau_{plm}$','$\tau_{tlm}$','$\tau_{tlc}$','Interpreter','latex','FontSize', FONTSIZE);
for idx=1:8
    icons(idx).LineWidth = 6;
end

grid on
hold off
%print('-dpng', '-r300', 'tau_master')
%%
figure 
hold on

plot(tau_pls(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(tau_tls(1:SUBSET), 'LineWidth', LINEWIDTH);
[~, icons] = legend('$\tau_{pls}$','$\tau_{tls}$','Interpreter','latex','FontSize', FONTSIZE);
for idx=1:6
    icons(idx).LineWidth = 6;
end

grid on
hold off

% [q_dot_ref,q_ref,qs_dot,qs,qm_dot,qm,qm_d,qs_d,qm_dot_d,qs_dot_d,fe_d,cnt,t,qs_d_prev,qm_d_prev,tau_plm,Hm,H-m,tau_tlc,tau_pls,Hs,H-s,fe,qm_m2s,qm_dot_m2s,tau_tls,qs_s2m,qs_dot_s2m,h_m_star,tau_tlm,H+m,H+s]
%%
figure
hold on
grid on

plot(Hm(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(Hs(1:SUBSET), 'LineWidth', LINEWIDTH);
[~, icons] = legend('$H_m$','$H_s$','Interpreter','latex','FontSize', FONTSIZE);
for idx=1:6
    icons(idx).LineWidth = 6;
end

hold off
%%
figure

hold on
grid on
plot(deltaH_s(1:SUBSET), 'LineWidth', LINEWIDTH);
plot(deltaH_m(1:SUBSET), 'LineWidth', LINEWIDTH);
[~, icons] = legend('$\Delta H_s$','$\Delta H_m$','Interpreter','latex','FontSize', FONTSIZE);
for idx=1:6
    icons(idx).LineWidth = 6;
end

hold off