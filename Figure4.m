clc
clear all

dynare NK
dynare NK_RR

model_RR = load('NK_RR_results')
model_TR = load('NK_results')

%% PC Slope

figure;

%subplot(2,1,1)
plot(model_RR.oo_.posterior_density.parameters.Kappa(:,1), model_RR.oo_.posterior_density.parameters.Kappa(:,2), '-k', 'linewidth', 4,'color',rgb('dimgray'))
hold on
xlabel('$\kappa$','interpreter','LaTeX','fontsize',22);
xlim([0 0.1])
ylim([0 inf])
%title('Real Rate Rule Model','interpreter','LaTeX','fontsize',22);
print('PCSlope2','-dpng');

figure;

%subplot(2,1,2)
plot(model_TR.oo_.posterior_density.parameters.Kappa(:,1), model_TR.oo_.posterior_density.parameters.Kappa(:,2), '-k', 'linewidth', 4,'color',rgb('dimgray'))
hold on
xlabel('$\kappa$','interpreter','LaTeX','fontsize',22);
xlim([0 1])
ylim([0 inf])
%title('Taylor Rule Model','interpreter','LaTeX','fontsize',22);

print('PCSlope1','-dpng');
