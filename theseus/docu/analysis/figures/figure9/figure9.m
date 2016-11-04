

topp=0.95;
w=0.33;
h=0.15;
tr = logspace(-2,8,50);



xspot=0.05;

subplot('Position',[xspot topp-h w h])
tc = 100000;
ta = 10;
nr = 1./(1+(tc./tr).^(3/4));
sr = 1./(tr/ta + ta./tr);

loglog(tr,sr,':','DisplayName','Stress')
hold on
loglog(tr,nr/2,':','DisplayName','Viscosity')
ylim([0 1])
xlim([10^-2 10^8])
xlabel('Recycling Time (\tau_r)')
set(gca,'ytick',[],'xtick',[ta tc], 'xticklabel',{'\tau_a' '\tau_c'},'xminortick','off')

semilogx(tr,sr./nr/10000,'Color',[0.25 0.25 0.25],'DisplayName','Strain Rate')

ylabel('Relative Units')
xlabel('Recycling Time (\tau_r)')
set(gca,'ytick',[],'xtick',[ta tc], 'xticklabel',{'\tau_a' '\tau_c'},'xminortick','off')
legend('Location','southwest')




subplot('Position',[xspot topp-h*2-0.05 w h])

tc = 10;
ta = 100000;
nr = 1./(1+(tc./tr).^(3/4));
sr = 1./(tr/ta + ta./tr);

loglog(tr,sr,':','DisplayName','Stress')
hold on
loglog(tr,nr/2,':','DisplayName','Viscosity')
ylim([0 1])
xlim([10^-2 10^8])

semilogx(tr,sr./nr/10,'Color',[0.25 0.25 0.25],'DisplayName','Strain Rate')

ylabel('Relative Units')
xlabel('Recycling Time (\tau_r)')
set(gca,'ytick',[],'xtick',[tc ta], 'xticklabel',{'\tau_c' '\tau_a'},'xminortick','off')






annotation('textbox', [0.01 0.9+0.025 0.05 0.05],'String','a)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.01 0.9-h-0.025 0.05 0.05],'String','b)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])
annotation('textbox', [0.41 0.9 0.05 0.05],'String','c)','LineStyle','none','FontSize',16,'FontName','Times','Color',[0.25 0.25 0.25])


xspot=xspot+0.42;
subplot('Position',[xspot topp-h*2-0.025 0.8-w h*2])

r = logspace(-3,3,50);
[X,Y] = meshgrid(logspace(-3,3,50),logspace(-3,3,50));
% C = (1+X.^(-3/4))./(X./Y+Y./X);
C = (1+(Y./X).^(3/4))./(X+1./X);

contour(X,Y,log10(C),20)
hold on
plot(r,r./r,':','Color',[0.25 0.25 0.25])
plot(r./r,r,':','Color',[0.25 0.25 0.25])

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\tau_r/\tau_a')
ylabel('\tau_c/\tau_a')
h = colorbar;
ylabel(h, 'Strain Rate')
set(h,'YTick',[])
set(h,'fontsize',14);
ylim([0.001,1000])

cd('../../../figures')
print('-depsc','-r0',['figure_theor.eps']);
