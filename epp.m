#! /usr/bin/octave
#
# evaluate-partitioning.present.m
#
# Randal A. Koene, 20050324

LOADPATH=[LOADPATH,':',system('printf $HOME',1),'/octave/common-m//'];

load -ascii evaluate-partitioning.ascii

N = [netsizes netsizes2D netsizes3D];
M = [maxsegs maxsegs2D maxsegs3D];
LN = log2(N);
LM = log2(M);

gplot N using 1:2, N using 1:3, N using 1:4, N using 1:5
plot_wait('evaluate_partitioning_netsizes');
closeplot

gplot M using 1:2, M using 1:3, M using 1:4, M using 1:5
plot_wait('evaluate_partitioning_maxsegs');
closeplot

gplot LN using 1:2, LN using 1:3, LN using 1:4, LN using 1:5
plot_wait('evaluate_partitioning_lognetsizes');
closeplot

gplot LM using 1:2, LM using 1:3, LM using 1:4, LM using 1:5
plot_wait('evaluate_partitioning_logmaxsegs');
closeplot

DN1 = diff(N(:,1));
DN2 = diff(N(:,2));
DN3 = diff(N(:,3));
DN4 = diff(N(:,4));
DN5 = diff(N(:,5));

DN = netsizes(2:size(netsizes,1));

slope2Dbruteforce = DN3./DN1

average = mean(slope2Dbruteforce)

slope2Dpartitioned = DN2./DN1

average = mean(slope2Dpartitioned)

slope3Dbruteforce = DN5./DN1

average = mean(slope3Dbruteforce)

slope3Dpartitioned = DN4./DN1

average = mean(slope3Dpartitioned)

S = [DN slope2Dpartitioned slope2Dbruteforce slope3Dpartitioned slope3Dbruteforce];
gplot S using 1:2, S using 1:3, S using 1:4, S using 1:5
plot_wait('evaluate_partitioning_netsizes_slope');
closeplot

DLN1 = diff(LN(:,1));
DLN2 = diff(LN(:,2));
DLN3 = diff(LN(:,3));
DLN4 = diff(LN(:,4));
DLN5 = diff(LN(:,5));

DLN = LN(2:size(netsizes,1),1);

slopelog2Dbruteforce = DLN3./DLN1

average = mean(slopelog2Dbruteforce)

slopelog2Dpartitioned = DLN2./DLN1

average = mean(slopelog2Dpartitioned)

slopelog3Dbruteforce = DLN5./DLN1

average = mean(slopelog3Dbruteforce)

slopelog3Dpartitioned = DLN4./DLN1

average = mean(slopelog3Dpartitioned)

S = [DLN slopelog2Dpartitioned slopelog2Dbruteforce slopelog3Dpartitioned slopelog3Dbruteforce];
gplot S using 1:2, S using 1:3, S using 1:4, S using 1:5
plot_wait('evaluate_partitioning_lognetsizes_slope');
closeplot
