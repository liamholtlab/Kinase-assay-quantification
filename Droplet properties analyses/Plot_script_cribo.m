%Figure plot script with errorbars under different conditions
load('Ribosome_concentration_calculation_data_new.mat')

y = pElk_hypershift_tubulin_rela;
dy = dpElk_hypershift_tubulin_rela;
% y = [pElk_hypershift_tubulin_rela_dhog1,pElk_hypershift_tubulin_rela_dhog1dsfp1];
% dy = [dpElk_hypershift_tubulin_rela_dhog1,dpElk_hypershift_tubulin_rela_dhog1dsfp1];
% y = [x_mCherry(1:4)./x_mCherry(1),x_mCherry(5:8)./x_mCherry(5)];
% dy = [dx_mCherry(1:4)./x_mCherry(1),dx_mCherry(5:8)./x_mCherry(5)];
% y = [x_GFP(1:4)./x_GFP(1),x_GFP(5:8)./x_GFP(5)];
% dy = [dx_GFP(1:4)./x_GFP(1),dx_GFP(5:8)./x_GFP(5)];
% y = [size_ratio(1:4)/size_ratio(1),size_ratio(5:8)/size_ratio(5)];
% dy = [dsize_ratio(1:4)/size_ratio(1),dsize_ratio(5:8)/size_ratio(5)];

x = c_ribo_rela * c_ribo0;
dx = dc_ribo_rela * c_ribo0;

N = length(x);
h = zeros(1,N);
C = {[0.00,0.45,0.74],[0.85,0.33,0.10],[0.93,0.69,0.13],[0.49,0.18,0.56],[0.64,0.08,0.18],[1.00,0.41,0.16],[0.30,0.75,0.93],[0.72,0.27,1.00]}; %Color array

figure
hold on
for i = 1:N
    h_temp = plot(x(i),y(i),'*','MarkerSize',8,'color',C{i});
    errorbar(x(i),y(i),dx(i),'horizontal','color',C{i})
    errorbar(x(i),y(i),dy(i),'color',C{i})
    
    h(i) = h_temp;
end
legend(h,name,'Location','northwest')

set(gca,'FontSize',16)

ylabel('Relative hypershift/tubulin')
% ylabel('Relative mean GFP-ERK1 intensity')
% ylabel('Relative mean mCherry-ELK1 intensity')
% ylabel('Relative droplet area fold change')

% xlabel('Relative ribosome concentration')
xlabel('Ribosome concentration (\muM)')