load('build/bin/values.mat')

figure
hold on

plot(qm_dot_ref)
plot(qs_dot)
plot(qm_dot)

grid on
title('Reference vs actual velocity trajectory')
hold off

// print('-dpng', '-r300', 'vel')

figure
hold on

plot(qm_ref)
plot(qs)
plot(qm)

grid on
title('Reference vs actual position trajectory')
hold off

// print('-dpng', '-r300', 'pos')