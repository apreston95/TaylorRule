clc 
clear all
close all

dynare SW2007_AP_RR_Full_Rule_IRFs noclearall;
dynare SW2007_IRFs noclearall;

%% 

H = length(IRF_OUTPUT_TFP_FULL_RR)-1;

figure;

subplot(7,3,1)
plot(0:H,IRF_OUTPUT_TFP_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_OUTPUT_TFP_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('TFP','interpreter','LaTeX','fontsize',8);
title('Output','interpreter','LaTeX','fontsize',12);
legend('SR','TR','interpreter','LaTeX','fontsize',8,'location','south','orientation','horizontal');
legend('boxoff')

subplot(7,3,2)
plot(0:H,IRF_INFLATION_TFP_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_INFLATION_TFP_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
title('Inflation','interpreter','LaTeX','fontsize',12);

subplot(7,3,3)
plot(0:H,IRF_NOMINAL_TFP_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_NOMINAL_TFP_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
yticks([-0.04 0 0.04])
ylim([-0.05 0.05])
line([0 H],[0 0],'color','black','LineStyle','-');
title('Nominal Interest Rate','interpreter','LaTeX','fontsize',12);






subplot(7,3,4)
plot(0:H,IRF_OUTPUT_RP_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_OUTPUT_RP_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Risk Premium','interpreter','LaTeX','fontsize',8);

subplot(7,3,5)
plot(0:H,IRF_INFLATION_RP_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_INFLATION_RP_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');

subplot(7,3,6)
plot(0:H,IRF_NOMINAL_RP_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_NOMINAL_RP_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');


subplot(7,3,7)
plot(0:H,IRF_OUTPUT_G_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_OUTPUT_G_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Government','interpreter','LaTeX','fontsize',8);

subplot(7,3,8)
plot(0:H,IRF_INFLATION_G_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_INFLATION_G_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');

subplot(7,3,9)
plot(0:H,IRF_NOMINAL_G_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_NOMINAL_G_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');



subplot(7,3,10)
plot(0:H,IRF_OUTPUT_I_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_OUTPUT_I_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Investment ','interpreter','LaTeX','fontsize',8);

subplot(7,3,11)
plot(0:H,IRF_INFLATION_I_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_INFLATION_I_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');


subplot(7,3,12)
plot(0:H,IRF_NOMINAL_I_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_NOMINAL_I_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');


subplot(7,3,13)
plot(0:H,IRF_OUTPUT_M_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_OUTPUT_M_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Monetary','interpreter','LaTeX','fontsize',8);

subplot(7,3,14)
plot(0:H,IRF_INFLATION_M_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_INFLATION_M_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');

subplot(7,3,15)
plot(0:H,IRF_NOMINAL_M_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_NOMINAL_M_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');



subplot(7,3,16)
plot(0:H,IRF_OUTPUT_PINF_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_OUTPUT_PINF_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Price Markup','interpreter','LaTeX','fontsize',8);

subplot(7,3,17)
plot(0:H,IRF_INFLATION_PINF_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_INFLATION_PINF_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');

subplot(7,3,18)
plot(0:H,IRF_NOMINAL_PINF_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_NOMINAL_PINF_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');


subplot(7,3,19)
plot(0:H,IRF_OUTPUT_W_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_OUTPUT_W_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Wage Markup','interpreter','LaTeX','fontsize',8);
xlabel('Quarters','interpreter','LaTeX','fontsize',8);

subplot(7,3,20)
plot(0:H,IRF_INFLATION_W_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_INFLATION_W_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
ylim([-0.1 0.15])
yticks([-0.1 0 0.1])
line([0 H],[0 0],'color','black','LineStyle','-');
xlabel('Quarters','interpreter','LaTeX','fontsize',8);

subplot(7,3,21)
plot(0:H,IRF_NOMINAL_W_FULL_RR(1:end),'color',rgb('dimgray'),'LineWidth',2);
hold on
plot(0:H,IRF_NOMINAL_W_TR(1:end),'color',rgb('lightgray'),'LineWidth',2);
hold on
axis auto
line([0 H],[0 0],'color','black','LineStyle','-');
yticks([0 0.03 0.06])
xlabel('Quarters','interpreter','LaTeX','fontsize',8);

print('IRFs_SW_RR','-dpng');

