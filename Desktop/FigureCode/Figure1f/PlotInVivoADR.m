%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plots for the in vivo ADR data from the bicistronic GCaMP-P2A-mRuby3
%Joey Broussard
%Tian Lab, UC Davis
%03/15/2017
%
%As of 02/10/2027, the data is structured with the ""ADR variable
%containing the normalized ADR values for all analyzed cells. Axon and soma
%have the green (first) and red (second) values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open the data
load('/FigureCode/Figure1f/ADRdata.mat')

%Plot the data
%get medians
G6sADRmed = median(G6sADR,'omitnan');
GAPADRmed = median(GAPADR,'omitnan');
%make the plot
figure;
histogram(G6sADR,'Normalization','pdf','Binwidth',0.2);
hold on;
histogram(GAPADR,'Normalization','pdf','Binwidth',0.22);
%add median lines and values
ylims = ylim;
plot([G6sADRmed,G6sADRmed],ylims,'b');
G6sADRmedText = num2str(G6sADRmed);
GAPADRmedText = num2str(GAPADRmed);
text(G6sADRmed,4*(ylims(2)-ylims(1))/5,G6sADRmedText,'color','b');
plot([GAPADRmed,GAPADRmed],ylims,'r');
text(GAPADRmed,3*(ylims(2)-ylims(1))/5,GAPADRmedText,'color','r');
xlabel('nADR')
ylabel('PDF')
legend('u-GCaMP6','axon-GCaMP6')
box off;
%adjust fontsize
set(findall(gcf,'-property','Fontsize'),'Fontsize',18);