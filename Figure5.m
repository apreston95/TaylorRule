clear all
close all

dynare NK_RR_PreVolcker
dynare NK_RR_PostVolcker

%% 


load('Determinacy_Cond_Post.mat')
load('Determinacy_Cond_Pre.mat')

figure;
plot(xi_Post,f_Post,'color',rgb('dimgray'),'LineWidth',4);
% xline(mean_Post,'--r');
hold on
plot(xi_Pre,f_Pre,'color',rgb('lightgray'),'LineWidth',4);
xline(1,'--r','color','black','LineWidth',2);
% xline(mean_Pre,'--r');
xlabel('$\Omega$','interpreter','LaTeX','fontsize',18);
ylabel('Posterior Density','interpreter','LaTeX','fontsize',18);
xlim([0 4]);
title('Determinacy Condition Posterior Distribution','interpreter','LaTeX','fontsize',18);
legend('Post-Volcker','Pre-Volcker','interpreter','LaTeX','fontsize',15);

print('Determinacy_Condition','-dpng');

%% 
figure;
ecdf(Determinacy_Post)

pts = 0:0.001:max(Determinacy_Pre);

[F_Post,Xi_Post] = ksdensity(Determinacy_Post,pts,'Function','cdf');
[F_Pre,Xi_Pre] = ksdensity(Determinacy_Pre,pts,'Function','cdf');

figure;

hold on
plot(Xi_Post,F_Post,'color',[0, 0, 1],'LineWidth',2)
hold on
plot(Xi_Pre,F_Pre,'color',[1, 0.5, 0],'LineWidth',2)
legend('Post-Volcker','Pre-Volcker','interpreter','LaTeX','fontsize',12);
xlabel('Determinacy Condition')
ylabel('Estimated cdf')

