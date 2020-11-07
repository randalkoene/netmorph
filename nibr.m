#! /usr/bin/octave
# Randal A. Koene, 20050415
# nibr.m - nibr statistics

clear

LOADPATH=[LOADPATH,':',system('printf $HOME',1),'/octave/common-m//'];

# Statistical data of simulated network development

if (version=="2.1.57")
  for i=1:size(argv,1),
    LOADPATH=[LOADPATH,':',nth(argv,i)];
  endfor
else
  if (nargin>0)
  for i=1:size(argv,1),
    LOADPATH=[LOADPATH,':',argv(i,:)];
  endfor
  endif
endif

stats

# In addition to statistical data, the following variables are
# now defined: filebase, autoplot

# For backward compatibility:
statsoutputcols = size(Dlength,2);
if (statsoutputcols>4)
  DlengthN = Dlength(:,1); Dlength = Dlength(:,2:statsoutputcols);
  AlengthN = Alength(:,1); Alength = Alength(:,2:statsoutputcols);
  DtermsegsperarborN = Dtermsegsperarbor(:,1); Dtermsegsperarbor = Dtermsegsperarbor(:,2:statsoutputcols);
  AtermsegsperarborN = Atermsegsperarbor(:,1); Atermsegsperarbor = Atermsegsperarbor(:,2:statsoutputcols);
  DterminallengthN = Dterminallength(:,1); Dterminallength = Dterminallength(:,2:statsoutputcols);
  AterminallengthN = Aterminallength(:,1); Aterminallength = Aterminallength(:,2:statsoutputcols);
  DintermediatelengthN = Dintermediatelength(:,1); Dintermediatelength = Dintermediatelength(:,2:statsoutputcols);
  AintermediatelengthN = Aintermediatelength(:,1); Aintermediatelength = Aintermediatelength(:,2:statsoutputcols);
  DlengthbetweenbifurcationsN = Dlengthbetweenbifurcations(:,1); Dlengthbetweenbifurcations = Dlengthbetweenbifurcations(:,2:statsoutputcols);
  AlengthbetweenbifurcationsN = Alengthbetweenbifurcations(:,1); Alengthbetweenbifurcations = Alengthbetweenbifurcations(:,2:statsoutputcols);
  DtermlensincebifurcationN = Dtermlensincebifurcation(:,1); Dtermlensincebifurcation = Dtermlensincebifurcation(:,2:statsoutputcols);
  AtermlensincebifurcationN = Atermlensincebifurcation(:,1); Atermlensincebifurcation = Atermlensincebifurcation(:,2:statsoutputcols);
  DturnsbetweenbifurcationsN = Dturnsbetweenbifurcations(:,1); Dturnsbetweenbifurcations = Dturnsbetweenbifurcations(:,2:statsoutputcols);
  AturnsbetweenbifurcationsN = Aturnsbetweenbifurcations(:,1); Aturnsbetweenbifurcations = Aturnsbetweenbifurcations(:,2:statsoutputcols);
  DbranchanglesN = Dbranchangles(:,1); Dbranchangles = Dbranchangles(:,2:statsoutputcols);
  AbranchanglesN = Abranchangles(:,1); Abranchangles = Abranchangles(:,2:statsoutputcols);
  DturnanglesN = Dturnangles(:,1); Dturnangles = Dturnangles(:,2:statsoutputcols);
  AturnanglesN = Aturnangles(:,1); Aturnangles = Aturnangles(:,2:statsoutputcols);
endif

# Computing remaining statistics of the simulation

meanidx=1; stdidx=2; minidx=3; maxidx=4;
ASQRTN = sqrt(Aarborssampled);
DSQRTN = sqrt(Darborssampled);
Dtwotailed_t = t_inv((1.0-(0.05/2.0))*ones(size(Darborssampled,1),1),Darborssampled-1);
Atwotailed_t = t_inv((1.0-(0.05/2.0))*ones(size(Aarborssampled,1),1),Aarborssampled-1);
global DNT = [T Dtermsegsperarbor];
DNT_ERR=Dtermsegsperarbor(:,stdidx) ./ DSQRTN;
DNT_CFL=Dtermsegsperarbor(:,meanidx) - (Dtwotailed_t.*DNT_ERR);
DNT_CFU=Dtermsegsperarbor(:,meanidx) + (Dtwotailed_t.*DNT_ERR);
DNT_ERR= [T DNT_ERR];
DNT_CFL= [T DNT_CFL];
DNT_CFU= [T DNT_CFU];
global ANT = [T Atermsegsperarbor];
ANT_ERR=Atermsegsperarbor(:,stdidx) ./ ASQRTN;
ANT_CFL=Atermsegsperarbor(:,meanidx) - (Atwotailed_t.*ANT_ERR);
ANT_CFU=Atermsegsperarbor(:,meanidx) + (Atwotailed_t.*ANT_ERR);
ANT_ERR= [T ANT_ERR];
ANT_CFL= [T ANT_CFL];
ANT_CFU= [T ANT_CFU];
global DL = [T Dlength];
DL_ERR=Dlength(:,stdidx) ./ DSQRTN;
DL_CFL=Dlength(:,meanidx) - (Dtwotailed_t.*DL_ERR);
DL_CFU=Dlength(:,meanidx) + (Dtwotailed_t.*DL_ERR);
DL_ERR= [T DL_ERR];
DL_CFL= [T DL_CFL];
DL_CFU= [T DL_CFU];
global AL = [T Alength];
AL_ERR=Alength(:,stdidx) ./ ASQRTN;
AL_CFL=Alength(:,meanidx) - (Atwotailed_t.*AL_ERR);
AL_CFU=Alength(:,meanidx) + (Atwotailed_t.*AL_ERR);
AL_ERR= [T AL_ERR];
AL_CFL= [T AL_CFL];
AL_CFU= [T AL_CFU];
global DTL = [T Dterminallength];
global ATL = [T Aterminallength];
global DIL = [T Dintermediatelength];
global AIL = [T Aintermediatelength];
global DLBB= [T Dlengthbetweenbifurcations];
global ALBB= [T Alengthbetweenbifurcations];
global DTLB= [T Dtermlensincebifurcation];
global ATLB= [T Atermlensincebifurcation];
global DTBB= [T Dturnsbetweenbifurcations];
global ATBB= [T Aturnsbetweenbifurcations];
global ABA = [T Abranchangles];
global DBA = [T Dbranchangles];
global ATA = [T Aturnangles];
global DTA = [T Dturnangles];
global A7BA = [T Asevenmicronbranchangles];
global D7BA = [T Dsevenmicronbranchangles];

# Experimental data

global total_axon_length;
global branch_points_per_axon;
global terminal_segments_per_axon;
global axons_mean_terminal_segment_length;
global axons_mean_intermediate_segment_length;
global number_of_dendrites_per_neuron;
global total_dendrite_length;
global branch_points_per_dendrite;
global terminal_segments_per_dendrite;
global dendrites_mean_terminal_segment_length;
global dendrites_mean_intermediate_segment_length;

printf("Loading culture experiment data for comparative measures\n");

load -ascii culture-experiment-data.ascii

global Eageidx=1; global Enidx=2; global Emeanidx=3; global Estdidx=4; global Esterridx=5; global Ecfiloidx=6; global Ecfihiidx=7; global Eminidx=8; global Emaxidx=9;
TE = total_axon_length(:,Eageidx);

printf("  [+ Loading branch angle sample data]\n");

load -ascii ramakers-branching-samples.ascii

# ID codes are of the form XXXXX[0/1], where XXXXX identifies a neuron and
# a terminating 0 indicates that this is data from a dendrite, while a
# terminating 1 indicates that this is data from an axon.

# Branch ID (BRID)
# Thickness of main element before branch (TMBB)
# Thickness of main element after branch (TMAB)
# Thickness of side element (TSE)
# Angle from main element to main element (Analpha)
# Angle from side element to main element (Anbeta)
# Radial distance to soma (RDTS)

global BRIDidx=1; global TMBBidx=2; global TMABidx=3; global TSEidx=4; global Analphaidx=5; Anbetaidx=6; RDTSidx=7;

# Here, convert the (currently single developmental age) data to statistical
# data points

global dendrites_seven_micron_branch_angles;
global axons_seven_micron_branch_angles;
global E7ageidx = 1;
global E7meanidx = 2;

printf("  [+ Loading turn angle sample data]\n");

load -ascii ramakers-turning-samples.ascii

# Compare simulation and experiment data visually

printf("Plotting network statistics\n");

if (autoplot==1)
     standard_auto_plots();
else

gplot terminal_segments_per_dendrite using Eageidx:Emeanidx, DNT using 1:2
plot_wait([filebase,'-dendritent']);
closeplot
gplot terminal_segments_per_axon using Eageidx:Emeanidx, ANT using 1:2
plot_wait([filebase,'-axonnt']);
closeplot
gplot total_dendrite_length using Eageidx:Emeanidx, DL using 1:2
plot_wait([filebase,'-dendritelength']);
closeplot
gplot total_axon_length using Eageidx:Emeanidx, AL using 1:2
plot_wait([filebase,'-axonlength']);
closeplot

gplot terminal_segments_per_dendrite using Eageidx:Estdidx, DNT using 1:3
plot_wait([filebase,'-dendritentstd']);
closeplot
gplot terminal_segments_per_axon using Eageidx:Estdidx, ANT using 1:3
plot_wait([filebase,'-axonntstd']);
closeplot
gplot total_dendrite_length using Eageidx:Estdidx, DL using 1:3
plot_wait([filebase,'-dendritelengthstd']);
closeplot
gplot total_axon_length using Eageidx:Estdidx, AL using 1:3
plot_wait([filebase,'-axonlengthstd']);
closeplot

gplot terminal_segments_per_dendrite using Eageidx:Emeanidx, DNT using 1:2, terminal_segments_per_dendrite using Eageidx:Ecfiloidx, terminal_segments_per_dendrite using Eageidx:Ecfihiidx, DNT_CFL using 1:2, DNT_CFU using 1:2
plot_wait([filebase,'-dendritentcfi']);
closeplot
gplot terminal_segments_per_axon using Eageidx:Emeanidx, ANT using 1:2, terminal_segments_per_axon using Eageidx:Ecfiloidx, terminal_segments_per_axon using Eageidx:Ecfihiidx, ANT_CFL using 1:2, ANT_CFU using 1:2
plot_wait([filebase,'-axonntcfi']);
closeplot
gplot total_dendrite_length using Eageidx:Emeanidx, DL using 1:2, total_dendrite_length using Eageidx:Ecfiloidx, total_dendrite_length using Eageidx:Ecfihiidx, DL_CFL using 1:2, DL_CFU using 1:2
plot_wait([filebase,'-dendritelengthcfi']);
closeplot
gplot total_axon_length using Eageidx:Emeanidx, AL using 1:2, total_axon_length using Eageidx:Ecfiloidx, total_axon_length using Eageidx:Ecfihiidx, AL_CFL using 1:2, AL_CFU using 1:2
plot_wait([filebase,'-axonlengthcfi']);
closeplot

# Test fitted function (this is just a quick output sample created by
# computing a fit with the lsqply program)

tstart = floor(T(1));
tend = floor(T(size(T,1)));
ALfittest = zeros(tend,3);
DLfittest = zeros(tend,3);
ANTfittest = zeros(tend,3);
DNTfittest = zeros(tend,3);
for t = tstart:tend,
  AlensmoothedO1 = -3227.891988 + (1210.271817*t);
  #AlensmoothedO2 = -1189.333527 + (577.490539*t) - 30.051510*(t*t);
  AlensmoothedO3 = -1846.00292 + (943.503247*t) - 13.712198*(t*t) + 1.343053*(t*t*t);
  ALfittest(t,1) = t;
  ALfittest(t,2) = AlensmoothedO1;
  #ALfittest(t,3) = AlensmoothedO2;
  ALfittest(t,3) = AlensmoothedO3;
  DlensmoothedO1 = -5.606921 + 6.089600*t;
  DlensmoothedO3 = 7.785732 + (3.538558*t) - 0.140114*(t*t) + 0.013301*(t*t*t);
  DLfittest(t,1) = t;
  DLfittest(t,2) = DlensmoothedO1;
  DLfittest(t,3) = DlensmoothedO3;
  AntsmoothedO1 = 3.127941 + 4.844813*t;
  AntsmoothedO3 = 3.455951 + (6.119460*t) - 0.284442*(t*t) + 0.011399*(t*t*t);
  ANTfittest(t,1) = t;
  ANTfittest(t,2) = AntsmoothedO1;
  ANTfittest(t,3) = AntsmoothedO3;
  DntsmoothedO1 = 1.135301 + 0.145156*t;
  DntsmoothedO3 = 1.422402 + (0.056846*t) + 0.004062*(t*t) + 0.000007*(t*t*t);
  DNTfittest(t,1) = t;
  DNTfittest(t,2) = DntsmoothedO1;
  DNTfittest(t,3) = DntsmoothedO3;
endfor

# Test fitted with mirroring around initial values (origin for length, 0-x and 2-y for number of terminal segments)
# (See TL#200505031005.2 for values used with lsqply.)

ALmirrored = zeros(tend,3);
DLmirrored = zeros(tend,3);
ANTmirrored = zeros(tend,3);
DNTmirrored = zeros(tend,3);
for t = tstart:tend,
  AlensmoothedO1 = 0.0 + (946.700573*t); # standard deviation computed by lsqply 3166.066301
  AlensmoothedO3 = 0.0 + (581.080071*t) - 0.0 + (1.339061*(t*t*t)); # standard deviation computed by lsqply 2226.362798
  ALmirrored(t,1) = t;
  ALmirrored(t,2) = AlensmoothedO1;
  ALmirrored(t,3) = AlensmoothedO3;
  DlensmoothedO1 = 0.0 + (5.631771*t); # standard deviation computed by lsqply 12.793726
  DlensmoothedO3 = 0.0 + (3.839628*t) - 0.0 + (0.006564*(t*t*t)); # standard deviation computed by lsqply 6.475500
  DLmirrored(t,1) = t;
  DLmirrored(t,2) = DlensmoothedO1;
  DLmirrored(t,3) = DlensmoothedO3;
  ANTsmoothedO1 = 1.0 + (5.018568*t); # standard deviation computed by lsqply 10.154034
  ANTsmoothedO3 = 1.0 + (4.724645*t) - 0.0 + (0.001076*(t*t*t)); # standard deviation computed by lsqply 9.991477
  ANTmirrored(t,1) = t;
  ANTmirrored(t,2) = ANTsmoothedO1;
  ANTmirrored(t,3) = ANTsmoothedO3;
  DNTsmoothedO1 = 1.0 + (0.156204*t); # standard deviation computed by lsqply 0.564477
  DNTsmoothedO3 = 1.0 + (0.145511*t) - 0.0 + (0.000039*(t*t*t)); # standard deviation computed by lsqply 0.560625
  DNTmirrored(t,1) = t;
  DNTmirrored(t,2) = DNTsmoothedO1;
  DNTmirrored(t,3) = DNTsmoothedO3;
endfor

gplot total_axon_length using Eageidx:Emeanidx, ALfittest using 1:2, ALfittest using 1:3, ALmirrored using 1:2, ALmirrored using 1:3, AL using 1:2
plot_wait([filebase,'-axonlengthsmoothedO1O3']);
closeplot
gplot total_dendrite_length using Eageidx:Emeanidx, DLfittest using 1:2, DLfittest using 1:3, DLmirrored using 1:2, DLmirrored using 1:3, DL using 1:2
plot_wait([filebase,'-dendritelengthsmoothedO1O3']);
closeplot
gplot terminal_segments_per_axon using Eageidx:Emeanidx, ANTfittest using 1:2, ANTfittest using 1:3, ANTmirrored using 1:2, ANTmirrored using 1:3, ANT using 1:2
plot_wait([filebase,'-axonntmoothedO1O3']);
closeplot
gplot terminal_segments_per_dendrite using Eageidx:Emeanidx, DNTfittest using 1:2, DNTfittest using 1:3, DNTmirrored using 1:2, DNTmirrored using 1:3, DNT using 1:2
plot_wait([filebase,'-dendritentmoothedO1O3']);
closeplot

endif

# Extract corresponding ages in the simulation to calculate fit

ALTE = AL(TE,:);
DLTE = DL(TE,:);

