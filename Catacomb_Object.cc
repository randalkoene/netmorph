/*
  © Copyright 2008 Randal A. Koene <randalk@netmorph.org>
  
  With design assistance from J. van Pelt & A. van Ooyen, and support
  from the Netherlands Organization for Scientific Research (NWO)
  Program Computational Life Sciences grant CLS2003 (635.100.005) and
  from the EC Marie Curie Research and Training Network (RTN)
  NEURoVERS-it 019247.

  This file is part of NETMORPH.

  NETMORPH is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  NETMORPH is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with NETMORPH.  If not, see <http://www.gnu.org/licenses/>.
*/
// Catacomb_Object.cc
// Randal A. Koene, 20080316

#include "diagnostic.hh"
#include "global.hh"
#include "network.hh"
#include "spatial.hh"
#include "Catacomb_Object.hh"

Catacomb_Group * ccmregionsgroup = NULL;
Catacomb_Group * ccmprojectionsgroup = NULL;
Catacomb_Group * ccmconnectionsgroup = NULL;
long Catacomb_neuronindex = 0;
long Catacomb_connectionindex = 0;

Catacomb_Header::Catacomb_Header(): text("<org.enorg.catacomb.core.CcmbBox>\n\
<ccmb_version>\"2.210\"</ccmb_version>\n\
<objects>\n\
</objects>\n\
<lists>\n\
<org.enorg.catacomb.core.CcmbList>\n\
<className>\"org.enorg.catacomb.lab.Workbench\"</className>\n\
<name>\"workbench\"</name>\n\
<org.enorg.catacomb.lab.Workbench>\n\
<assemblyItems>\n\
<className>\"org.enorg.catacomb.assembly.AssemblyItem\"</className>\n\
<name>\"assemblyItems\"</name>\n\
<org.enorg.catacomb.run.FixedStepRunner>\n\
<runtime>100.0</runtime>\n\
<timestep>0.2</timestep>\n\
<maxTcalcByTreal>0.1</maxTcalcByTreal>\n\
<updateFrequency>10.0</updateFrequency>\n\
<autoRerun>0</autoRerun>\n\
<position>\n\
  -17.27623406341942  8.591042824452437 \n\
</position>\n\
<scale>\n\
  1.0  1.0 \n\
</scale>\n\
<color>10065101</color>\n\
<locked>0</locked>\n\
<name>\"fixedStepRunner\"</name>\n\
</org.enorg.catacomb.run.FixedStepRunner>\n") {
}

String Catacomb_Neuron::str() {
  // Note that due to the normal ordering used for output,
  // the output for this object begins with the closing
  // of the preceding workbench content.
  String s("</assemblyItems>\n\
<color>14465878</color>\n\
<locked>0</locked>\n\
<name>\"Workbench\"</name>\n\
</org.enorg.catacomb.lab.Workbench>\n\
</org.enorg.catacomb.core.CcmbList>\n\
<org.enorg.catacomb.core.CcmbList>\n\
<className>\"org.enorg.catacomb.network.Neuron\"</className>\n\
<name>\"neuron\"</name>\n\
<org.enorg.catacomb.network.Neuron>\n\
<optimizerHint>1</optimizerHint>\n\
<assemblyItems>\n\
<className>\"org.enorg.catacomb.assembly.AssemblyItem\"</className>\n\
<name>\"assemblyItems\"</name>\n\
<org.enorg.catacomb.cell.IntegratorCompartment>\n\
<radius>5.0</radius>\n\
<membraneCapacitance>1.0</membraneCapacitance>\n\
<leak>2</leak>\n\
<resistance>1.0</resistance>\n\
<timeConstant>5.0</timeConstant>\n\
<spikeWidth>1.0</spikeWidth>\n\
<resetPotential>-70.0</resetPotential>\n\
<restPotential>-70.0</restPotential>\n\
<position>\n\
-1.1237332461588743  -1.2827721477607072 \n\
</position>\n\
<scale>\n\
1.0  1.0 \n\
</scale>\n\
<color>14197955</color>\n\
<locked>0</locked>\n\
<name>\"integratorCompartment\"</name>\n\
</org.enorg.catacomb.cell.IntegratorCompartment>\n\
<org.enorg.catacomb.cell.SpikeResponseFunction>\n\
<Erev>-90.0</Erev>\n\
<conductance>1.0</conductance>\n\
<riseTime>4.0</riseTime>\n\
<fallTime>8.0</fallTime>\n\
<duration>30.0</duration>\n\
<position>\n\
0.006047728015692755  1.4004576659038896 \n\
</position>\n\
<scale>\n\
1.0  1.0 \n\
</scale>\n\
<color>5741236</color>\n\
<locked>0</locked>\n\
<name>\"spikeResponseFunction\"</name>\n\
</org.enorg.catacomb.cell.SpikeResponseFunction>\n\
<org.enorg.catacomb.cell.ThresholdSpikeGenerator>\n\
<threshold>-50.0</threshold>\n\
<refractoryPeriod>1.0</refractoryPeriod>\n\
<position>\n\
2.3362209872507362  -0.5766590389016022 \n\
</position>\n\
<scale>\n\
1.0  1.0 \n\
</scale>\n\
<color>8857063</color>\n\
<locked>0</locked>\n\
<name>\"thresholdSpikeGenerator\"</name>\n\
</org.enorg.catacomb.cell.ThresholdSpikeGenerator>\n\
<org.enorg.catacomb.cell.SynapsePopulation>\n\
<Erev>0.0</Erev>\n\
<G_0>0.1</G_0>\n\
<G_max>10.0</G_max>\n\
<riseTime>4.0</riseTime>\n\
<fallTime>10.0</fallTime>\n\
<duration>30.0</duration>\n\
<saturates>0</saturates>\n\
<saturationFactor>0.5</saturationFactor>\n\
<voltageDependent>0</voltageDependent>\n\
<Vhalf>0.0</Vhalf>\n\
<z>0.0</z>\n\
<position>\n\
-3.595129127165741  0.8355671788166061 \n\
</position>\n\
<scale>\n\
1.0  1.0 \n\
</scale>\n\
<color>16332687</color>\n\
<locked>0</locked>\n\
<name>\"synapsePopulation\"</name>\n\
</org.enorg.catacomb.cell.SynapsePopulation>\n\
<org.enorg.catacomb.link.InsertionArrow>\n\
<position>\n\
-2.4888852566198096  -0.5531219352729657 \n\
</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n\
\"-0.3662	0.6486	synapsePopulation	membrane	-1\"\n\
\"0.7151	-0.07965	integratorCompartment	channels	-1\"\n\
</portData>\n\
<color>12162413</color>\n\
<locked>0</locked>\n\
<name>\"insertionArrow\"</name>\n\
</org.enorg.catacomb.link.InsertionArrow>\n\
<org.enorg.catacomb.link.InsertionArrow>\n\
<position>\n\
-0.27639751552795033  0.15299117358613845 \n\
</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n\
\"0.2824	0.2474	spikeResponseFunction	membrane	-1\"\n\
\"-1.497	-0.7857	integratorCompartment	channels	-1\"\n\
</portData>\n\
<color>9685697</color>\n\
<locked>0</locked>\n\
<name>\"insertionArrow_0\"</name>\n\
</org.enorg.catacomb.link.InsertionArrow>\n\
<org.enorg.catacomb.link.VectorCable>\n\
<position>\n\
0.7592350441320708  -0.6472703497875125 \n\
</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n\
\"-1.383	-0.1355	integratorCompartment	potential	-1\"\n\
\"0.5769	0.07061	thresholdSpikeGenerator	potentialIn	-1\"\n\
</portData>\n\
<color>7065045</color>\n\
<locked>0</locked>\n\
<name>\"vectorCable\"</name>\n\
</org.enorg.catacomb.link.VectorCable>\n\
<org.enorg.catacomb.link.SpikeCable>\n\
<position>\n\
2.595129127165741  2.106570774762994 \n\
</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n\
\"0.7410	-2.683	thresholdSpikeGenerator	spikeOut	-1\"\n\
\"-2.589	0.2938	spikeResponseFunction	spikeIn	-1\"\n\
</portData>\n\
<color>15739752</color>\n\
<locked>0</locked>\n\
<name>\"spikeCable\"</name>\n\
</org.enorg.catacomb.link.SpikeCable>\n\
<org.enorg.catacomb.link.SpikeCable>\n\
<position>\n\
2.38329519450801  -2.106570774762995 \n\
</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n\
\"0.9529	1.529	thresholdSpikeGenerator	spikeOut	-1\"\n\
\"-2.857	0.1738	integratorCompartment	reset	-1\"\n\
</portData>\n\
<color>13846473</color>\n\
<locked>0</locked>\n\
<name>\"spikeCable_0\"</name>\n\
</org.enorg.catacomb.link.SpikeCable>\n\
<org.enorg.catacomb.cell.SynapsePopulation>\n\
<Erev>0.0</Erev>\n\
<G_0>0.1</G_0>\n\
<G_max>10.0</G_max>\n\
<riseTime>4.0</riseTime>\n\
<fallTime>10.0</fallTime>\n\
<duration>30.0</duration>\n\
<saturates>0</saturates>\n\
<saturationFactor>0.5</saturationFactor>\n\
<voltageDependent>0</voltageDependent>\n\
<Vhalf>0.0</Vhalf>\n\
<z>0.0</z>\n\
<position>\n\
-3.586862886561444  -2.0404827899054783 \n\
</position>\n\
<scale>\n\
1.0  1.0 \n\
</scale>\n\
<color>6223864</color>\n\
<locked>0</locked>\n\
<name>\"synapsePopulation_0\"</name>\n\
</org.enorg.catacomb.cell.SynapsePopulation>\n\
<org.enorg.catacomb.link.InsertionArrow>\n\
<position>\n\
-2.2376720897815927  -1.9099159386042022 \n\
</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n\
\"-0.6091	-0.8705	synapsePopulation_0	membrane	-1\"\n\
\"0.4639	1.277	integratorCompartment	channels	-1\"\n\
</portData>\n\
<color>12823593</color>\n\
<locked>0</locked>\n\
<name>\"insertionArrow_1\"</name>\n\
</org.enorg.catacomb.link.InsertionArrow>\n\
</assemblyItems>\n\
<color>11372167</color>\n\
<locked>0</locked>\n\
<name>\"neuron\"</name>\n\
</org.enorg.catacomb.network.Neuron>\n\
</org.enorg.catacomb.core.CcmbList>\n");
  return s;
}

String Catacomb_Region::str() {
  String s("<org.enorg.catacomb.lab.VectorRecorder>\n\
<sampleInterval>0.3</sampleInterval>\n\
<recordFromStart>1</recordFromStart>\n\
<events>1</events>\n\
<continuous>1</continuous>\n\
<plot>1</plot>\n\
<display>\n\
<shownItems>\n");
  for (int i=0; i<n; i++) s += "1  ";
  s += "\n</shownItems>\n\
<color>4911530</color>\n\
<locked>0</locked>\n\
<name>\"display\"</name>\n\
</display>\n\
<resultsInIcon>0</resultsInIcon>\n\
<position>\n";
  s += String(x+5.0,"%.5f  ") + String(y+5.0,"%.5f \n"); 
  s += "</position>\n\
<scale>\n\
1.0  1.0 \n\
</scale>\n\
<color>16678431</color>\n\
<locked>0</locked>\n\
<name>\"";
  s += name + "_recorder";
  s += "\"</name>\n\
</org.enorg.catacomb.lab.VectorRecorder>\n\
<org.enorg.catacomb.network.NeuronSocket>\n\
<neuron>\n\
<className>\"org.enorg.catacomb.network.Neuron\"</className>\n\
<selectedIndex>0</selectedIndex>\n\
<selectedName>\"neuron\"</selectedName>\n\
<name>\"neuron\"</name>\n\
</neuron>\n\
<singleElement>0</singleElement>\n\
<nrow>";
  s += String((long) n);
  s += "</nrow>\n\
<ncolumn>1</ncolumn>\n\
<independent>0</independent>\n\
<position>\n";
  s += String(x,"%.5f  ") + String(y,"%.5f \n");
  s += "</position>\n\
<scale>\n\
4.0  4.0 \n\
</scale>\n\
<color>2981602</color>\n\
<locked>0</locked>\n\
<name>\"";
  s += name;
  s += "\"</name>\n\
</org.enorg.catacomb.network.NeuronSocket>\n";
  s += "<org.enorg.catacomb.spike.JoinCloneRelay>\n\
<position>\n";
  s += String(x+4.0,"%.5f  ") + String(y+4.0,"%.5f \n"); 
  s += "</position>\n\
<scale>\n\
1.0  1.0 \n\
</scale>\n\
<color>15674009</color>\n\
<locked>0</locked>\n\
<name>\"";
  s += name + "_JoinCloneRelay";
  s += "\"</name>\n\
</org.enorg.catacomb.spike.JoinCloneRelay>\n";
  s += "<org.enorg.catacomb.link.SpikeCable>\n\
<position>\n";
  s += String(x+2.0,"%.5f  ") + String(y+2.0,"%.5f \n"); 
  s += "</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n";
  s += "\"-0.4416	-6.088	" + name + "	thresholdSpikeGenerator.spikeOut	-1\"\n\
\"0.5534	1.435	" + name + "_JoinCloneRelay	spikes in	0\"\n";
  s += "</portData>\n\
<color>10514627</color>\n\
<locked>0</locked>\n\
<name>\"";
  s += name + "toRelay_spikeCable\"</name>\n\
</org.enorg.catacomb.link.SpikeCable>\n";
  s += "<org.enorg.catacomb.link.SpikeCable>\n\
<position>\n";
  s += String(x+4.5,"%.5f  ") + String(y+4.5,"%.5f \n"); 
  s += "</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n";
  s += "\"-0.8829	-1.059	" + name + "_JoinCloneRelay	spike out	-1\"\n\
\"0.4122	0.7176	" + name + "_recorder	spike in	-1\"\n";
  s += "</portData>\n\
<color>1702589</color>\n\
<locked>0</locked>\n\
<name>\"";
  s += "Relayto" + name + "_recorder_spikeCable\"</name>\n\
</org.enorg.catacomb.link.SpikeCable>\n";
  s += "<org.enorg.catacomb.link.VectorCable>\n\
<position>\n";
  s += String(x+2.5,"%.5f  ") + String(y+2.5,"%.5f \n"); 
  s += "</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>16777215</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n";
  s += "\"-3.015	-4.364	" + name + "	thresholdSpikeGenerator.effective potential	-1\"\n\
\"1.189	2.977	" + name + "_recorder	vector in	-1\"\n";
  s += "</portData>\n\
<color>14930997</color>\n\
<locked>0</locked>\n\
<name>\"";
  s += name + "torecorder_vectorCable\"</name>\n\
</org.enorg.catacomb.link.VectorCable>\n";
  return s;
}

String Catacomb_Projection::str() {
  String s;
  for (int i=0; i<ccmprojections->Dim(); i++) {
    for (int j=0; j<ccmprojections->Dim(); j++) {
      if (ccmprojections->el(i,j)!=NULL) {
	s += "<org.enorg.catacomb.network.SpikeProjection>\n\
<mapping>3</mapping>\n\
<connectionTable>\n\
<className>\"org.enorg.catacomb.network.ConnectionTable\"</className>\n\
<selectedIndex>0</selectedIndex>\n\
<selectedName>\"";
	region * r_from = regions->el(i);
	region * r_to = regions->el(j);
	s += r_from->Name()+"to"+r_to->Name();
	s += "\"</selectedName>\n\
<name>\"connectionTable\"</name>\n\
</connectionTable>\n\
<applyTableWeights>0</applyTableWeights>\n\
<applyTableDelays>0</applyTableDelays>\n\
<applyDelay>0</applyDelay>\n\
<delay>0.0</delay>\n\
<position>\n";
	  spatial p(r_from->center());
	if (i==j) {
	  s += String(p.X()+5.0,"%.5f  ") + String(p.Y()-5.0,"%.5f\n");
	} else {
	  p += r_to->center();
	  p /= 2.0;
	  s += String(p.X(),"%.5f  ") + String(p.Y(),"%.5f\n");
	}
	s += "</position>\n\
<scale>\n\
2.0  2.0 \n\
</scale>\n\
<label>\"\"</label>\n\
<background>0</background>\n\
<foreground>0</foreground>\n\
<hideLabel>0</hideLabel>\n\
<portData>\n";
	s += "\"7.486	4.759	" + r_from->Name() + "	thresholdSpikeGenerator.spikeOut	-1\"\n";
	s += "\"0.3859	6.058	" + r_to->Name();
	if (i==j) s += "	synapsePopulation.afferent spike	-1\"\n";
	else s += "	synapsePopulation_0.afferent spike	-1\"\n";
	s += "</portData>\n\
<color>7124075</color>\n\
<locked>0</locked>\n\
<name>\"";
	s += r_from->Name()+"to"+r_to->Name() + "_spikeProjection\"</name>\n\
</org.enorg.catacomb.network.SpikeProjection>\n";
      }
    } 
  }
  return s;
}

String Catacomb_ConnectionTable::str() {
  String s("<org.enorg.catacomb.core.CcmbList>\n\
<className>\"org.enorg.catacomb.network.ConnectionTable\"</className>\n\
<name>\"connectionTable\"</name>\n");
  for (int i=0; i<ccmprojections->Dim(); i++) {
    for (int j=0; j<ccmprojections->Dim(); j++) {
      ccm_tableptr & ctp = ccmprojections->el(i,j); // a reference can be set only once
      if (ctp!=NULL) {
	s += "<org.enorg.catacomb.network.ConnectionTable>\n\
<selectable>1</selectable>\n\
<n_in>";
	region * r_from = regions->el(i);
	region * r_to = regions->el(j);
	s += String((long) ctp->Dim_in());
	s += "</n_in>\n\
<n_out>";
	s += String((long) ctp->Dim_out());
	s += "</n_out>\n\
<new_weight>0.8</new_weight>\n\
<new_delay>2.0</new_delay>\n\
<importer>\n\
<file>\"\"</file>\n\
<format>1</format>\n\
<color>13148976</color>\n\
<locked>0</locked>\n\
<name>\"importer\"</name>\n\
</importer>\n\
<table>\n\
<sourceIDs>\n\
</sourceIDs>\n\
<targetIDs>\n\
</targetIDs>\n\
<nrow>";
	s += String((long) ctp->Dim_in());
	s += "</nrow>\n\
<ncol>";
	s += String((long) ctp->Dim_out());
	s += "</ncol>\n\
<rows>\n";
	for (int r=0; r<ctp->Dim_in(); r++) {
	  s += "<org.enorg.sand.network.ConnectionRow>\n\
<sourceID>\"\"</sourceID>\n\
<icon>\n";
	  int numindices = 0;
	  for (int c=0; c<ctp->Dim_in(); c++) if (ctp->el(r,c)!=0) {
	    s += String((long) c);
	    numindices++;
	  }
	  s += "\n</icon>\n\
<weights>\n";
	  for (int ni=0; ni<numindices; ni++) s += "0.8 "; // [***NOTE] Change this when specifying weights from the NETMORPH abstract connection strengths.
	  s += "\n</weights>\n\
<delays>\n";
	  for (int ni=0; ni<numindices; ni++) s += "2.0 ";
	  s += "\n</delays>\n\
</org.enorg.sand.network.ConnectionRow>\n";
	}
	s += "</rows>\n\
</table>\n\
<color>6155905</color>\n\
<locked>0</locked>\n\
<name>\"";
	s += r_from->Name() + "to" + r_to->Name();
	s += "\"</name>\n\
</org.enorg.catacomb.network.ConnectionTable>\n";
      }
    }
  }
  s += "</org.enorg.catacomb.core.CcmbList>\n";
  return s;
}

String Catacomb_Group::str() {
  String s;
    if (!comment.empty()) {
      comment.prepend("<!-- ");
      comment += " -->";
      s = comment;
    }
  PLL_LOOP_FORWARD(Catacomb_Object,ccmobjects.head(),1) s+=e->str();
  return s;
}

void Catacomb_Group::Add_Catacomb_Object(Catacomb_Object * ccmobject) {
  if (!ccmobject) return;
  ccmobjects.link_before(ccmobject);
}

