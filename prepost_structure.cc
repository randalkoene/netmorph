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
// prepost_structure.cc
// Randal A. Koene, 20041118

#include "global.hh"
#include "prepost_structure.hh"
#include "axon_direction_model.hh"
#include "neuron.hh"
#include "Command_Line_Parameters.hh"

natural_schema_parent_set axons_most_specific_natural_set[UNTYPED_NEURON+1] = {
  all_axons_sps,
  all_interneuron_axons_sps,
  all_axons_sps,
  all_pyramidal_axons_sps,
  all_axons_sps
};

presynaptic_structure::presynaptic_structure(neuron & _n, Segment & initseg, spatial & acoords, int _typeid): fibre_structure(_n,initseg,acoords) {
  // If a general or specific parameter settings so indicate, the terminal
  // segment created in this constructor can be associated with a direction
  // model.
  // [***INCOMPLETE] Only general parameter settings are recognized so far.
  // See TL#200511060626.2, DIL#20051105072005.1 and TL#200603111137.1.
  // Despite this incompleteness, the linking of contributing models is already
  // done through arbor-specific detection in the CLP.
  // PROTOCOL: As described in task log entries of DIL#20070612153815.1, protocol-compliant sets have
  // model schema assigned to them. Therefore, components need only look for the schema of the most
  // specific (narrow) set to which they belong.
#ifndef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  axon_direction_history * adh = NULL;
  if (general_axon_direction_model!=legacy_dm) {
    switch (general_axon_direction_model) {
    case segment_history_tension_dm: 
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
      adh = new pll_axon_direction_history();
#else
      adh = new parse_axon_direction_history(*this);
#endif
      terminalsegments.head()->set_direction_model(new tension_direction_model(*adh,general_axon_direction_model_root,*main_clp));
      break;;
    case cell_attraction_dm:
      terminalsegments.head()->set_direction_model(new cell_attraction_direction_model(general_axon_direction_model_root,*main_clp));
      break;;
    default: break;; // to avoid a compilation warning message about legacy_dm
    }
  }
#endif // NOT(INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS)
  // [***INCOMPLETE] This is a crude way to parse applicable subsets
  // (see TL#200603120618.4)
  // Beware: I assume elmodel==NULL during initialization!
  // A hierarchy of natural subsets (see DIL#20060323042401.1)
  terminal_segment * ts = terminalsegments.head();
  neuron_type ntype = N()->TypeID();
  natural_schema_parent_set nset = axons_most_specific_natural_set[ntype];
  bmmodel = branching_model_subset_schemas[nset]->clone();
  ts->set_branch_angle_model(branch_angle_model_subset_schemas[nset]->clone());
  elmodel = general_arbor_elongation_model_base_subset_schemas[nset]->clone();
  ts->set_elongation_model(general_terminal_segment_elongation_model_base_subset_schemas[nset]->clone());
  ts->set_ERI_model(elongation_rate_initialization_model_subset_schemas[nset]->clone());
#ifdef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  ts->set_direction_model(direction_model_subset_schemas[nset]->clone());
#endif
}

terminal_segment_ptr * presynaptic_structure::array_of_terminal_segments() {
  // an array with pointers to all terminal segments in all axonal trees
  terminalarraylength = terminalsegments.length();
  PLL_LOOP_FORWARD(fibre_structure,Next(),1) terminalarraylength += e->TerminalSegments()->length();
  if (terminalarraylength==0) return NULL;
  terminal_segment_ptr * tsarray = new terminal_segment_ptr[terminalarraylength];
  int i = 0;
  PLL_LOOP_FORWARD_NESTED(fibre_structure,this,1,ps) {
    PLL_LOOP_FORWARD(terminal_segment,ps->TerminalSegments()->head(),1) tsarray[i++] = e;
  }
  return tsarray;
}

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
void presynaptic_structure::Collect_Data(network_statistics_base & nsb) {
  if (nsb.netstats[axons_termsegs_n].Collect()) nsb.netstats[axons_termsegs_n].Collect_Data((double) terminalsegments.length());
  fibre_collected_data fcd(nsb,true);
  PLL_LOOP_FORWARD(terminal_segment,terminalsegments.head(),1) e->Collect_Data(fcd); // this must be done first
  fibre_segment::Collect_Data(fcd);
  nsb.netstats[axons_arbor_len].Collect_Data(fcd.arborlength);
}
#endif

natural_schema_parent_set dendrites_most_specific_natural_set[UNTYPED_NEURON+1] = {
  all_dendrites_sps,
  all_interneuron_dendrites_sps,
  all_dendrites_sps,
  all_pyramidal_dendrites_sps,
  all_dendrites_sps
};

postsynaptic_structure::postsynaptic_structure(neuron & _n, Segment & initseg, spatial & acoords, int _typeid): fibre_structure(_n,initseg,acoords) {
  // If a general or specific parameter settings so indicate, the terminal
  // segment created in this constructor can be associated with a direction
  // model.
  // [***INCOMPLETE] Only general parameter settings are recognized so far.
  // See TL#200511060626.2.
  // PROTOCOL: As described in task log entries of DIL#20070612153815.1, protocol-compliant sets have
  // model schema assigned to them. Therefore, components need only look for the schema of the most
  // specific (narrow) set to which they belong.
#ifndef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  axon_direction_history * adh = NULL;
  if (general_dendrite_direction_model!=radial_dm) {
    switch (general_dendrite_direction_model) {
    case segment_history_tension_dm:
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
      adh = new pll_axon_direction_history();
#else
      adh = new parse_axon_direction_history(*this);
#endif
      terminalsegments.head()->set_direction_model(new tension_direction_model(*adh,general_dendrite_direction_model_root,*main_clp));
      break;;
   case cell_attraction_dm:
      terminalsegments.head()->set_direction_model(new cell_attraction_direction_model(general_dendrite_direction_model_root,*main_clp));
      break;
    default: break;; // to avoid a compilation warning message about legacy_dm
    }
  }
#endif
  // [***INCOMPLETE] This is a crude way to parse applicable subsets
  // (see TL#200603120618.4)
  // A hierarchy of natural subsets (see DIL#20060323042401.1)
  terminal_segment * ts = terminalsegments.head();
  neuron_type ntype = N()->TypeID();
  natural_schema_parent_set nset = dendrites_most_specific_natural_set[ntype];
  if ((ntype==PYRAMIDAL) && (_typeid==ALL_APICAL_PYRAMIDAL_DENDRITES_AEM)) nset = all_apical_pyramidal_dendrites_sps;
  bmmodel = branching_model_subset_schemas[nset]->clone();
  ts->set_branch_angle_model(branch_angle_model_subset_schemas[nset]->clone());
  elmodel = general_arbor_elongation_model_base_subset_schemas[nset]->clone();
  ts->set_elongation_model(general_terminal_segment_elongation_model_base_subset_schemas[nset]->clone());
  ts->set_ERI_model(elongation_rate_initialization_model_subset_schemas[nset]->clone());
#ifdef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  ts->set_direction_model(direction_model_subset_schemas[nset]->clone());
#endif

}

terminal_segment_ptr * postsynaptic_structure::array_of_terminal_segments() {
  // an array with pointers to all terminal segments in all dendritic trees
  terminalarraylength = terminalsegments.length();
  PLL_LOOP_FORWARD(fibre_structure,Next(),1) terminalarraylength += e->TerminalSegments()->length();
  if (terminalarraylength==0) return NULL;
  terminal_segment_ptr * tsarray = new terminal_segment_ptr[terminalarraylength];
  int i = 0;
  PLL_LOOP_FORWARD_NESTED(fibre_structure,this,1,ps) {
    PLL_LOOP_FORWARD(terminal_segment,ps->TerminalSegments()->head(),1) tsarray[i++] = e;
  }
  return tsarray;
}

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
void postsynaptic_structure::Collect_Data(network_statistics_base & nsb) {
  if (nsb.netstats[dendrites_termsegs_n].Collect()) nsb.netstats[dendrites_termsegs_n].Collect_Data((double) terminalsegments.length());
  fibre_collected_data fcd(nsb,false);
  PLL_LOOP_FORWARD(terminal_segment,terminalsegments.head(),1) e->Collect_Data(fcd); // this must be done first
  fibre_segment::Collect_Data(fcd);
  nsb.netstats[dendrites_arbor_len].Collect_Data(fcd.arborlength);
}
#endif

