% This script is used to plot rainbow tracks based on the tracks_.mat file
% with jet color map having the tracks start as blue and tracks end as red.

% Input the name of the tracks you want to plot.
clear
load('tracked_Traj_TS68(4h) hog1dsfp1d cytoGEM DMSO 1MSorbitol 004.nd2_crop.mat')

figure
hold on
for i = 1:size(result,1)
    x = result(i).tracking.x;
    y = result(i).tracking.y;
    c = (1:length(x))/length(x);
    if length(x)>2
        patch([x' NaN],[y' NaN],[c NaN],[c NaN],'edgecolor','interp','LineWidth',2)
    else
        plot(x,y,'b','LineWidth',2)
    end
end
colormap('jet')
set(gca,'YDir','reverse')
colorbar
    
