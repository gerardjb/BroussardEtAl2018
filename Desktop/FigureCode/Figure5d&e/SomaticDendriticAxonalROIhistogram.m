%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Making Histogram for ROIs from Green to Red ratio analysis
%Joey Broussard
%Tian Lab, UC Davis
%05/23/2017
%
%Plotting histograms of the ratios of green and red fluorescence in the
%various cellular compartments studied for the labeled Scnn1a-Tg3-Cre
%animals.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data
load('SomaticG2Rratio.mat')
load('DendriticG2Rratio.mat')
load('AxonG2Rratio.mat')

%Collect variables
%Somata
G6s = y_G6s19_J91_MAP2;
GAP = [y_GAP621_J88_MAP2;y_GAP621_J88_noStain];
%Axons
G6s_axon = [G6s0523_J88_MAP2_axon;G6s0619_J91_L1_Axon;G6sNoStain_J88_L1_Axon];
GAP_axon = [GAP0523_J91_Apical_axon;GAP0619_J88_AL1V1_axon];
%Dendrites
G6s_dendrite = [G6s0523_J88_MAP2_dendrite;G6sNoStain_J88_L4_dendrite;G6s0619_J91_L1_dendrite];
GAP_dendrite = [GAP0523_J91_L4_dendrite;GAP_NoStain_J91_L4_dendrite];

%Make the figures
%Somata
figure;
histogram(G6s,'Binwidth',0.05,'Normalization','PDF');hold on;
histogram(GAP,'Binwidth',0.05,'Normalization','PDF')
legend('G6s','GAP')
ylims = ylim;
plot([median(G6s),median(G6s)],[ylims(1),ylims(2)],'b')
text(median(G6s),3,num2str(median(G6s)),'Color','b')
plot([median(GAP),median(GAP)],[ylims(1),ylims(2)],'r')
text(median(GAP),3,num2str(median(GAP)),'Color','r')
title('Somatic green to red ratio')
xlabel('Somatic greent to red ratio (A.U.)')
ylabel('PDF');
%Axons
figure;
histogram(G6s_axon,'Binwidth',0.1,'Normalization','PDF');hold on;
histogram(GAP_axon,'Binwidth',0.1,'Normalization','PDF')
legend('G6s_axon','GAP_axon')
ylims = ylim;
plot([median(G6s_axon),median(G6s_axon)],[ylims(1),ylims(2)],'b')
text(median(G6s_axon),3,num2str(median(G6s_axon)),'Color','b')
plot([median(GAP_axon),median(GAP_axon)],[ylims(1),ylims(2)],'r')
text(median(GAP_axon),3,num2str(median(GAP_axon)),'Color','r')
title('Axonal green to red ratio')
xlabel('Axonal greent to red ratio (A.U.)')
ylabel('PDF');
%Dendrites
figure;
histogram(G6s_dendrite,'Binwidth',0.05,'Normalization','PDF');hold on;
histogram(GAP_dendrite,'Binwidth',0.05,'Normalization','PDF')
legend('G6s_dendrite','GAP_dendrite')
ylims = ylim;
plot([median(G6s_dendrite),median(G6s_dendrite)],[ylims(1),ylims(2)],'b')
text(median(G6s_dendrite),3,num2str(median(G6s_dendrite)),'Color','b')
plot([median(GAP_dendrite),median(GAP_dendrite)],[ylims(1),ylims(2)],'r')
text(median(GAP_dendrite),3,num2str(median(GAP_dendrite)),'Color','r')
title('Dendritic green to red ratio')
xlabel('Dendritic greent to red ratio (A.U.)')
ylabel('PDF');