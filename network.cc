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
// network.cc
// Randal A. Koene, 20041118

#include <math.h>
#include <iostream>
#include <cstdint>
#include <memory>
#include "global.hh"
#include "file.hh"
#include "Network_Generated_Statistics.hh"
#include "Sampled_Output.hh"
#include "synapse_formation_model.hh"
#include "neurite_diameter_model.hh"
#include "connection.hh"
#include "synapse_structure.hh"
#include "network.hh"
#include "neuron.hh"
#include "Include/Embeddable.h"

// Python embedded as a const char string
#include "objcurves.py"

// variables

bool numneurons_is_default = false; // flag to indicate that numneurons was set by a fallback to default

// classes

region & region::append(String _name, neuron & n) {
  if (_name==name) return append(n);
  if (head()) return head()->conditional_create(_name,n);
  else return seek_head()->conditional_create(_name,n);
}

region & region::conditional_create(String _name, neuron & n) {
  if (_name==name) return append(n);
  if (Next()) return Next()->conditional_create(_name,n);
  if (Root()) {
    region * reg = new region(_name,n);
    Root()->link_before(reg);
    return *reg;
  } else {
    warning("Error: Unable to add regions without a PLLRoot<region> object in region::conditional_create()!\n");
    return *this;
  } 
}

int region::find(neuron & n) {
  // Returns the element number of the neuron in a region.
  // Returns -1 if not found.
  int i=0;
  PLL_LOOP_FORWARD(neuronptrlist,nlist.head(),1) {
    if (e->N()==&n) return i;
    i++;
  }
  return -1;
}

spatial region::center() {
  spatial regioncenter;
  PLL_LOOP_FORWARD(neuronptrlist,nlist.head(),1) regioncenter += e->N()->Pos();
  regioncenter /= (double) nlist.length();
  return regioncenter;
}

Catacomb_Group * region::net_Catacomb() {
  spatial C(center());
#define CATACOMBSCALE 0.5
  ccmregionsgroup->Add_Catacomb_Object(new Catacomb_Region(name,nlist.length(),CATACOMBSCALE*C.X(),CATACOMBSCALE*C.Y()));
  return ccmregionsgroup;
}

int regionslist::find(neuron & n) {
  // Returns the element number of the region in a list of
  // regions, which contains the neuron.
  // Returns -1 if not found.
  int i=0;
  PLL_LOOP_FORWARD(region,head(),1) {
    if (e->find(n)>=0) return i;
    i++;
  }
  return -1;
}

int regionslist::find_by_name(const String& _name) {
  int i=0;
  PLL_LOOP_FORWARD(region,head(),1) {
    if (e->Name()==_name) return i;
    i++;
  }
  return -1;
}

network::network(int numneurons, bool e): Event_Queue(this), netinfo("initialized with "+String((long) numneurons)+" neurons\n"), edges(e), candidate_synapses(true), synapses_during_development(true), paramspreadmin(NULL), paramspreadmax(NULL), nsr(NULL), sss(NULL), ndm(NULL) {
  // creates untyped neurons
  netinfo += "principal neurons\n";
  if (edges) netinfo += "with functional edges (possible edge effects)\n";
  else netinfo += "no functional edges (connections loop around)\n";
  for (int i=0; i<numneurons; i++) PLLRoot<neuron>::link_before(new principal());
}

network::network(int numneurons, double principalprobability, bool e): Event_Queue(this), netinfo("initialized with "+String((long) numneurons)+" neurons\n"), edges(e), candidate_synapses(true), synapses_during_development(true), paramspreadmin(NULL), paramspreadmax(NULL), nsr(NULL), sss(NULL), ndm(NULL) {
  // creates principal neurons and interneurons
  if (principalprobability>1.0) principalprobability=1.0;
  netinfo += "random principal neurons and interneurons, with a probability "+String(principalprobability,"%.2f")+" to generate a principal neuron\n";
  if (edges) netinfo += "with functional edges (possible edge effects)\n";
  else netinfo += "no functional edges (connections loop around)\n";
  unsigned int threshold = (unsigned int) (((double) X_misc.get_max())*principalprobability);
  for (int i=0; i<numneurons; i++) {
    if (X_misc.get_rand()<=threshold) PLLRoot<neuron>::link_before(new principal());
    else PLLRoot<neuron>::link_before(new interneuron());
  }
}

network::network(int numneurons, neuron * psmin, neuron * psmax, double principalprobability, bool e): Event_Queue(this), netinfo("initialized with "+String((long) numneurons)+" neurons\n"), edges(e), candidate_synapses(true), synapses_during_development(true), paramspreadmin(psmin), paramspreadmax(psmax), nsr(NULL), sss(NULL), ndm(NULL) {
  // creates principal neurons and interneurons
  if (principalprobability>1.0) principalprobability=1.0;
  netinfo += "random principal neurona and interneurons, with a probability "+String(principalprobability,"%.2f")+" to generate a principal neuron\n";
  if (edges) netinfo += "with functional edges (possible edge effects)\n";
  else netinfo += "no functional edges (connections loop around)\n";
  unsigned int threshold = (unsigned int) (((double) X_misc.get_max())*principalprobability);
  bool varyradius = false; double radiusvariationfactor = 0.0;
  if ((psmin) && (psmax)) if (psmin->Radius()<psmax->Radius()) {
    varyradius=true;
    radiusvariationfactor = (psmax->Radius()-psmin->Radius())/((double) X_misc.get_max());
  }
  neuron * n; double r;
  for (int i=0; i<numneurons; i++) {
    if (X_misc.get_rand()<=threshold) n = new principal();
    else n = new interneuron();
    if (varyradius) {
      r = (radiusvariationfactor*((double) X_misc.get_rand())) + psmin->Radius();
	 n->set_radius(r);
    }
    PLLRoot<neuron>::link_before(n);
  }
}

String network::show_target_approach() {
  String res;
  for (auto & [np, npa] : neuron_pair_map) {
    res += "Approach from " + np.from->numerical_IDStr() + " to " + np.to->numerical_IDStr() + ": ";
    for (auto & sqd : npa.sqdhistory) {
      res += String(sqd,"%.3f ");
    }
    res += " last entry t=" + String(npa.t_sample,"%.1f\n");
  }
  return res;
}

void network::target_approach_samples(neuron * from, neuron * to, double sqdistance) {
  neuron_pair np(from, to);
  auto it = neuron_pair_map.find(np);
  if (it==neuron_pair_map.end()) {
    neuron_pair_map[np] = neuron_pair_approach(Time(), sqdistance);
    return;
  }
  if (sqdistance < it->second.sqdistance) {
    it->second.sqdistance = sqdistance;
    it->second.t_sample = Time();
    it->second.sqdhistory.emplace_back(sqdistance);
  }
}


void network::parse_CLP(Command_Line_Parameters & clp) {
  // This parses specifications of general parameters. Very specialized
  // parameters are parsed within relevant functions (see below).
  Event_Queue::parse_CLP(clp);
  int n;
  if ((n=clp.Specifies_Parameter("candidate_synapses"))>=0) candidate_synapses = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("synapses_during_development"))>=0) synapses_during_development = (downcase(clp.ParValue(n))==String("true"));
  ndm = select_neurite_diameter_model(*this,clp);

  if ((n=clp.Specifies_Parameter("NES_output"))>=0) NES_output = (downcase(clp.ParValue(n))==String("true"));

  if ((n=clp.Specifies_Parameter("track_approach"))>=0) track_approach = (downcase(clp.ParValue(n))==String("true"));

  if ((n=clp.Specifies_Parameter("detailed_chemical_factors"))>=0) chemdata.detailed_chemical_factors = (downcase(clp.ParValue(n))==String("true"));

  if ((n=clp.Specifies_Parameter("soma_attraction_weight"))>=0) chemdata.soma_attraction_weight = atof(clp.ParValue(n));

  if ((n=clp.Specifies_Parameter("enable_completion_requirements"))>=0) completion.enable_completion_requirements = (downcase(clp.ParValue(n))==String("true"));
}

// Find neurons specified by region+type+index. E.g. "Mid.pyramidal.2".
void network::_innerloop_byregiontypecell_neuron_specific_configurator(Command_Line_Parameters & clp, int i, String identifier, String valuestr) {
  enum idparts_enum: int {
    id_region,
    id_type,
    id_idx,
    id_request,
    NUM_idparts_enum
  };
  StringList idparts(identifier,'.');
  if (idparts.length()<int(NUM_idparts_enum)) {
    warning("Missing region / cell type / index / request at command line parameter: "+clp.ParName(i));
    return;
  }

  if ((idparts.length()>int(NUM_idparts_enum)) && (!idparts[NUM_idparts_enum].empty())) {
    warning("Neuron identifier parsing with multipart requests not yet implemented: "+clp.ParName(i));
    return;
  }

  long idx = atol(idparts[id_idx].chars());

  PLL_LOOP_FORWARD(region, Regions().head(), 1) {
    if (downcase(e->Name())==idparts[id_region]) { // Correct region.
        
      long celltype_idx = 0;
      PLL_LOOP_FORWARD_NESTED(neuronptrlist, e->Nlist().head(), 1, _nptr) {
        neuron* n = _nptr->N();

        if (n->TypeStr()==idparts[id_type]) { // Correct cell type.

          if (celltype_idx==idx) { // Correct index.

            n->neuron_specific_configurator(idparts[id_request], valuestr);
            return;
          }
          celltype_idx++;
        }
      }

      return;  
    }
  }

  warning("Match not found for: "+clp.ParName(i));

}

void network::_innerloop_byidx_neuron_specific_configurator(Command_Line_Parameters & clp, int i, String identifier, String valuestr) {
  String neuron_idxstr(identifier.before('.'));
  String config_clp(identifier.after('.'));
  if (config_clp.empty()) {
    warning("Missing configuration request at command line parameter: "+clp.ParName(i));
  } else {
    neuron * n = find_neuron_by_idx(atol(neuron_idxstr.chars()));
    n->neuron_specific_configurator(config_clp, clp.ParValue(i));
  }
}

void network::_innerloop_neuron_specific_configurator(Command_Line_Parameters & clp, int i, String identifier, String valuestr) {
  if (identifier.empty()) {
    warning("Missing identifier details after 'idx.' in "+clp.ParName(i)+'\n');
    return;
  }

  // Find neurons specified by index in list of all neurons.
  char first = identifier[0];
  if ((first>='0') && (first<='9')) {
    _innerloop_byidx_neuron_specific_configurator(clp, i, identifier, valuestr);
  } else {

    _innerloop_byregiontypecell_neuron_specific_configurator(clp, i, identifier, valuestr);
  }

}

/**
 * Check for explicit "input=" declaration to identify input receiving
 * neurons. These are use for I/O path specific region-to-region connectome
 * completion checks. If not specified explicitly, these are later
 * determined automatically by `deduce_input_neurons_by_chemattraction()`.
 */
void network::explicit_input_neurons(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("inputs"))>=0) {
    String inputsstr(clp.URI_unescape_ParValue(n));
    while (!inputsstr.empty()) {
      // Find a region name
      String regionname(inputsstr.before(' '));
      if (regionname.empty()) {
        regionname = inputsstr;
        inputsstr = "";
      } else {
        inputsstr = inputsstr.after(' ');
      }

      // Find all the neurons in that region
      int num_region = Region_by_name(regionname);
      if (num_region<0) {
        warning("Unknown region: "+regionname);
      } else {
        region* r = Region_by_idx(num_region);
        PLL_LOOP_FORWARD(neuronptrlist, r->Nlist().head(), 1) {
          input_neurons.emplace_back(e->N());
        }
      }
    }
  }
  if (input_neurons.size()<1) error("No input neurons specified.\n");
  else {
    report("Input neurons: \n");
    for (auto& nptr : input_neurons) {
      report(nptr->label()+'\n');
    }
  }
}

/**
 * This applies a neuron identification protocol within an already-known network
 * architecture to identify specific neurons and to call a configuration function
 * on those.
 * 
 * Identification examples are:
 * "idx.0.attractors=" - specify chemical attraction factors emitted by neuron with index 0.
 * "idx.Mid.pyramidal.2.attractedto=" - specify chemical factors that the third pyramidal
 * neuron in region Mid is attracted to.
 */
void network::neuron_specific_configurator(Command_Line_Parameters & clp) {
  for (int i=0; i<clp.NumParameters(); i++) {
    // Find all the arguments that begin with "idx.".
    if (downcase(clp.ParName(i).before(4))=="idx.") {

      _innerloop_neuron_specific_configurator(clp, i, downcase(clp.ParName(i).after(3)), clp.ParValue(i));
      
    }
  }
}

void completion_requirements::update_connection(neuron* presyn, neuron* postsyn) {
  neuron_pair np(presyn, postsyn);
  auto it = target_connections.find(np);
  if (it != target_connections.end()) {
    it->second = true;
    //progress("Completion requirements: "+String(completion_percentage(), "%.1f\n"));
  }
}

String completion_requirements::report_completion_criteria() {
  String res("Target Connections:\n");
  for (auto& [np, state] : target_connections) {
    res += np.from->label() + '(' + np.from->numerical_IDStr() + ") -> " + np.to->label() + '(' + np.to->numerical_IDStr() + ")\n";
  }
  return res;
}

String completion_requirements::report_criteria_completed(network* net) {
  String res("Achieved Target Connections:\n");
  for (auto& [np, state] : target_connections) {
    if (state) {
      res += np.from->label() + '(' + np.from->numerical_IDStr() + ") -> " + np.to->label() + '(' + np.to->numerical_IDStr() + ")\n";
    }
  }
  if (net->input_neurons.empty()) net->input_neurons = net->deduce_input_neurons_by_chemattraction();
  double percentage = net->compare_region_connection_requirements(net->input_neurons);
  res += "Achieved Target Region Connections (" + String(percentage, "%.1f") + "%):\n";
  for (auto& [rp, state] : region_target_connections) {
    if (state>0) {
      res += net->Region_by_idx(rp.from)->Name() +" (" + String(long(state))+ ") -> " + net->Region_by_idx(rp.to)->Name() + '\n';
    }
  }
  return res;
}

void completion_requirements::criteria_numbers(long& totnum, long& completed) {
  totnum = 0;
  completed = 0;
  // Check target connections:
  totnum += target_connections.size();
  for (auto& [np, iscomplete] : target_connections) {
    if (iscomplete) completed++;
  }
}

double completion_requirements::completion_percentage() {
  long totnum, completed;
  criteria_numbers(totnum, completed);
  return 100.0*double(completed)/double(totnum);
}

bool completion_requirements::is_complete() {
  if (!enable_completion_requirements) return false;
  long totnum, completed;
  criteria_numbers(totnum, completed);
  return completed >= totnum;
}

void completion_requirements::region_criteria_numbers(long& totnum, long& completed) {
  totnum = 0;
  completed = 0;
  // Check target connections:
  totnum += region_target_connections.size();
  for (auto& [rp, num_connections] : region_target_connections) {
    if (num_connections>0) completed++;
  }
}

double completion_requirements::region_completion_percentage() {
  long totnum, completed;
  region_criteria_numbers(totnum, completed);
  return 100.0*double(completed)/double(totnum);
}

/**
 * This prepares the completion structure.
 * E.g. pairs of neurons that are attractors and attractedto are
 * used to set up target_connections.
 * Also, pairs of regions that are attractors and attractedto are
 * used to set up region_target_connections.
 * 
 * *** NOTE: In a more advanced version, completion requirements
 *     could be that at least one source neuron is connected to
 *     a target neuron, or something like that, or some minimum
 *     number.
 */
void network::setup_completion_requirements() {
  // Prepare the target_connections map.
  PLL_LOOP_FORWARD(neuron, PLLRoot<neuron>::head(), 1) {
    if (!e->chemdata.attractedto.empty()) {
      int region_from = regions.find(*e);
      for (auto& chemfactor : e->chemdata.attractedto) {
        if (chemdata.attractor_somata.find(chemfactor)==chemdata.attractor_somata.end()) {
          error("Missing attractor for factor "+chemindex_to_label[chemfactor]);
        } else {
          for (auto& target_nptr : chemdata.attractor_somata.at(chemfactor)) {
            neuron_pair np(e, target_nptr);
            completion.target_connections[np] = false;
            int region_to = regions.find(*const_cast<neuron*>(target_nptr));
            if ((region_from>=0) && (region_to>=0)) {
              region_number_pair rp(region_from, region_to);
              completion.region_target_connections[rp] = 0;
            }
          }
        }
      }
    }
  }
}

class has_no_chemattraction: public neuron_list_op {
protected:
  std::vector<neuron*>& nptr_set;
public:
  has_no_chemattraction(std::vector<neuron*>& _nptr_set): nptr_set(_nptr_set) {}
  virtual void op(neuron* n) {
    if (!n->is_attractor()) {
      nptr_set.emplace_back(n);
    }
  }
};

/**
 * Returns a set of neurons identified as input neurons based on
 * having no specified attractor chemicals.
 */
std::vector<neuron*> network::deduce_input_neurons_by_chemattraction() {
  std::vector<neuron*> input_neurons;
  has_no_chemattraction nonchemattractors(input_neurons);
  neuron_op(nonchemattractors);
  return input_neurons;
}

/**
 * Use this to set an integer cache value in neurons.
 */
class neuron_set_cache_int: public neuron_list_op {
public:
  int value = 0;
  neuron_set_cache_int() {}
  neuron_set_cache_int(int _value): value(_value) {}
  virtual void op(neuron* n) {
    n->cache.i = value;
  }
};

/**
 * From a non-traversed neuron traverse all connections, find target
 * neurons and traverse further.
 * While traversing, check all region pairs of the neurons involved
 * and if they are in the region_target_connections list then set
 * that pair's completion to true.
 * The main activity for this op all takes place in the
 * connection_op() call.
 */
class region_connections_search: public synapse_tree_op {
protected:
  network* net = nullptr;
public:
  region_connections_search(network* _net): net(_net) {}
  // This is run whenever entering a new neuron, if started there or above.
  virtual void neuron_op(neuron* n) {
    n->cache.i = 1; // Set the traversed flag for this neuron.
  }
  // This is run whenever entering a new connection, if started there or above.
  virtual void connection_op(connection * c) {
    // At the neuron level, this is running through the OutputConnections, from
    // pre to post.
    neuron* post_traverse = c->PostSynaptic();
    int region_from = net->Regions().find(*c->PreSynaptic());
    int region_to = net->Regions().find(*c->PostSynaptic());
    if ((region_from>=0) && (region_to>=0)) {
      region_number_pair rp(region_from, region_to);
      auto it = net->completion.region_target_connections.find(rp);
      if (it != net->completion.region_target_connections.end()) {
        it->second++;
      }
    }
    if (post_traverse->cache.i == 0) post_traverse->synapse_op(*this);
  }
  // This is run whenever entering a new synapse.
  virtual void op(synapse* s) {};
};

/**
 * Using a specified set of input neurons, search recursively
 * through connected neurons to compare established region-to-region
 * connection pairs with the set of region-to-region connections
 * required in region_target_connections.
 * Indirect, unspecified connection paths are normally permitted,
 * but only connection paths that begin at the input neurons are
 * checked, because it is assumed that possible I/O functions are
 * the reason for connection requirements.
 * Returns the percentage of required region connections completed.
 */
double network::compare_region_connection_requirements(std::vector<neuron*> input_neurons) {
  // Prepare cache flags of all neurons in the network.
  neuron_set_cache_int clear_flags(0);
  neuron_op(clear_flags);
  // Carry out this search for all of the neurons in the specified
  // input_regions.
  region_connections_search regconnsearch(this);
  for (auto& nptr : input_neurons) {
    nptr->synapse_op(regconnsearch);
  }
  return completion.region_completion_percentage();
}

String network::report_parameters() {
  String res(Event_Queue::report_parameters());
  res += "General network parameters:\n";
  if (candidate_synapses) res += "  Seek candidate synapse sites.\n";
  else res += "  Do not seek candidate synapse sites.\n";
  if (synapses_during_development) res += "  Generate actual synapses from candidates DURING development.\n";
  else res += "  Generate actual synapses from candidates AFTER development.\n";
  if (ndm) res += "  " + ndm->report_parameters() + '\n';
  return res;
}

String network::neuron_specific_reports() {
  String res;
  PLL_LOOP_FORWARD(neuron, PLLRoot<neuron>::head(), 1) {
    res += e->neuron_specific_reports();
  }
  return res;
}

void network::add_typed_neuron(neuron_type nt, neuron * psmin, neuron * psmax) {
  // adds a specific type of neuron to the network and prepares its
  // parameters
  principal * n;
  switch (nt) {
  case PRINCIPAL_NEURON: n = new principal(); break;
  case INTERNEURON: n = new interneuron(); break;
  case MULTIPOLAR_NONPYRAMIDAL: n = new multipolar_nonpyramidal(); break;
  case BIPOLAR: n = new bipolar(); break;
  case PYRAMIDAL: n = new pyramidal(); break;
  default: n = new principal(); break; //neuron(); break;
  }
  n->parse_CLP(*main_clp);
  if ((psmin) && (psmax)) if (psmin->Radius()<psmax->Radius()) {
     n->set_radius(X_misc.get_rand_range_real1(psmin->Radius(),psmax->Radius()));
   }
  // [***NOTE] Here I insure that neurons are randomly distributed, even when
  // exact population sizes are specified.
  int len = PLLRoot<neuron>::length();
  if (len<1) PLLRoot<neuron>::link_before(n);
  else {
    int pos = X_misc.get_rand_range(0,len);
    if (pos>=len) PLLRoot<neuron>::link_before(n);
    else PLLRoot<neuron>::insert_before(pos,n);
  }
}

void network::remove_abstract_connections_without_synapses() {
  // This function can be used to clean up abstract connections in networks
  // that require actual synapses to be defined for connections.
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    connection * c_next;
    for (connection * c = e->OutputConnections()->head(); (c); c = c_next) {
      c_next = c->Next(); // cached, since c may be removed
      if (!c->Synapses()->head()) c->remove();
    }
  }
}

network::network(int numneurons, neuron * psmin, neuron * psmax, Network_Statistics_Root & netstats, bool _edges): Event_Queue(this), netinfo("initialized with "+String((long) numneurons)+" general population neurons\n"), edges(_edges), candidate_synapses(true), synapses_during_development(true), paramspreadmin(psmin), paramspreadmax(psmax), nsr(&netstats), sss(NULL), ndm(NULL) {
  // creates multiple types of neurons
  netinfo += "multiple neuron types with the following network statistics where\nratios are values in a probability distribution:\n";
  netinfo += netstats.head()->all_str() + '\n';
  if (edges) netinfo += "with functional edges (possible edge effects)\n";
  else netinfo += "no functional edges (connections loop around)\n";
  int numneurontypes = netstats.length();
  if (numneurontypes<=0) error("nibr Error: Empty Network_Statistics in network::network()\n");
  if (netstats.UseExactPopulationSizes()) {
    PLL_LOOP_FORWARD(Network_Statistics,netstats.head(),1) for (int i = e->PopulationSize(); i>0; i--) add_typed_neuron(e->TypeID(),psmin,psmax);
  } else {
    // Alternatively, allocate by proportions if there are any general population neurons
    if (numneurons<1) return;
    int * cumprobdist = new int[numneurontypes];
    neuron_type * ntypesarray = new neuron_type[numneurontypes];
    int netstatidx = 0, cumprob = 0; double sumratio = 0.0;
    PLL_LOOP_FORWARD(Network_Statistics,netstats.head(),1) sumratio += e->Ratio();
    PLL_LOOP_FORWARD(Network_Statistics,netstats.head(),1) {
      ntypesarray[netstatidx] = e->TypeID();
      cumprob += (int) ((e->Ratio()/sumratio)*((double) X_misc.get_max()));
      cumprobdist[netstatidx] = cumprob;
      netstatidx++;
    }
    cumprobdist[numneurontypes-1] = X_misc.get_max();
    int nt_r, nt_idx;
    for (int i=0; i<numneurons; i++) {
      nt_r = X_misc.get_rand();
      for (nt_idx=0; nt_idx<numneurontypes; nt_idx++) if (cumprobdist[nt_idx]>=nt_r) break;
      add_typed_neuron(ntypesarray[nt_idx],psmin,psmax);
    }
    delete[] cumprobdist; delete[] ntypesarray;
  }
}

int network::get_or_add_chemlabel(String chemlabel) {
  auto it = chemlabel_to_index.find(chemlabel);
  if (it == chemlabel_to_index.end()) {
    int new_idx = chemindex_to_label.size();
    chemlabel_to_index[chemlabel] = new_idx;
    chemindex_to_label.emplace_back(chemlabel);
    return new_idx;
  } else {
    return it->second;
  }
}

/**
 * If this fiber segment is dendritic, and if the associated cell is an
 * attractor, then the P1 location of the new_branching_segment becomes
 * a new attractor point.
 * 
 * *** Note: At present, this assigns the same chemical factors as those
 *     used by the soma.
 */
void network::possibly_add_attractor(fibre_segment * new_branching_segment) {
  if (chemdata.has_specified_factors && chemdata.detailed_chemical_factors) {
    if (new_branching_segment->get_fstype()==dendrite_fs) {
      neuron * nptr = new_branching_segment->N();
      if (nptr->is_attractor()) {
        for (auto & chemfactor : nptr->chemdata.attractor) {
          chemdata.attractor_segments[chemfactor].emplace(new_branching_segment);
        }
      }
    }
  }
}

Shape_Hexagon_Result network::shape_hexagon(Command_Line_Parameters & clp) {
// gives neurons in the network positions such that they form a hexagon and
// are evenly distributed.
// sidelength is the length of each side of the hexagon in micrometers
  Shape_Hexagon_Result res;
  double sidelength = 10.0;
  int n;
  if ((n=clp.Specifies_Parameter("shape_horizontal"))>=0) sidelength = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("shape_sidelength"))>=0) sidelength = atof(clp.ParValue(n));
  res.displaywidth = 2.0*sidelength;
  int numneurons = PLLRoot<neuron>::length();
  if (numneurons<1) return res;
  netinfo += "shape hexagon with "+String(sidelength,"%.2f")+" micrometer sides\n";
  // hegaxon = 6 equilateral triangles
  // rectangular area = 1.5*length of side of triangle * 2*height of triangle
  // 1. distribute evenly in rectangular area so that all edges are also
  // populated by neurons
  // 2. if a neuron has an x position beyond the lines defining the right hand
  // boundary of the hexagon then give it a new x position boundary_x - x
  double triangleheight = sqrt((sidelength*sidelength) - (0.25*sidelength*sidelength));
  res.rectangulararea = (1.5*sidelength) * (2*triangleheight);
  netinfo += "area "+String(res.rectangulararea,"%.3f")+" micrometers^2\n";
  double areaperneuron = res.rectangulararea/((double) numneurons);
  double squaresidelength = sqrt(areaperneuron);
  res.neuronsperline = (int) (1.5*sidelength/squaresidelength);
  if (res.neuronsperline<1) res.neuronsperline = 1;
  // NOTE: After the (int) approximation above, preceding calculations are
  // approximations! For instance, we do not use squaresidelength directly
  // for the neuron spacing below for that reason.
  int discardneurons = numneurons % res.neuronsperline;
  if (discardneurons>0) {
    for (int i = 0; i<discardneurons; i++) PLLRoot<neuron>::tail()->remove();
    numneurons = PLLRoot<neuron>::length();
    res.message = "Number of neurons reduced for even distribution.\nNew number of neurons = " + String((long) numneurons) + "\n";
  }
  netinfo += "shaped network consists of "+String((long) numneurons)+" neurons\n";
  double rectangleneuronxspacing = (1.5*sidelength)/((double) (res.neuronsperline-1)); // see reason in NOTE above
  double rectangleneuronyspacing = (2.0*triangleheight)/((double) ((numneurons/res.neuronsperline)-1)); // see reason in NOTE above
  netinfo += "horizontal distance between neurons is "+String(rectangleneuronxspacing,"%.3f")+" micrometers\n";
  netinfo += "vertical distance between neurons is "+String(rectangleneuronyspacing,"%.3f")+" micrometers\n";
  double lineoneslope = sidelength/(2.0*triangleheight);
  double lineoneoffset = sidelength;
  double linetwoslope = -lineoneslope;
  double linetwooffset = 2.0*sidelength;
  int colnum = 0; double liney = 0.0;
  double boundaryonex = lineoneoffset, boundarytwox = linetwooffset, boundaryx = linetwooffset;
  if (lineoneoffset<linetwooffset) boundaryx=lineoneoffset;
  double lineonediff = lineoneslope*rectangleneuronxspacing;
  double linetwodiff = linetwoslope*rectangleneuronxspacing;
  double neuronx, neurony; double colval = 0.0;
  double halfhexagonneuronsoffset = rectangleneuronyspacing * ((double) (((numneurons/res.neuronsperline)-1)/2));
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    if (colval<0.0) { // mirroring part of rectangle back to create a hexagon
      colval -= 1.0;
      if (liney<halfhexagonneuronsoffset) neurony = liney+halfhexagonneuronsoffset;
      else neurony = liney - halfhexagonneuronsoffset;
      neuronx = colval*rectangleneuronxspacing;
    } else { // creating rectangle with the same area
      neurony = liney;
      neuronx = ((double) colnum)*rectangleneuronxspacing;
      if (neuronx>boundaryx) {
	colval = -1.0;
	neuronx = -rectangleneuronxspacing;
	if (liney<halfhexagonneuronsoffset) neurony = liney+halfhexagonneuronsoffset;
	else neurony = liney - halfhexagonneuronsoffset;
      }
    }
    e->set_position_in_Z_plane(neuronx,neurony);
    colnum++;
    if (colnum>=res.neuronsperline) {
      colnum=0; colval = 0.0;
      liney += rectangleneuronyspacing;
      boundaryonex += lineonediff;
      boundarytwox += linetwodiff;
      if (boundaryonex<boundarytwox) boundaryx=boundaryonex;
      else boundaryx=boundarytwox;
    }
  }
  // 3. insure that no more than one neuron occupies one position
  // Allocate to one network region
  neuron * rootneuron = PLLRoot<neuron>::head();
  region & reg = regions.append("net",*rootneuron);
  if (rootneuron) PLL_LOOP_FORWARD(neuron,rootneuron->Next(),1) reg.append(*e);
  res.message += "Number of neurons per line = "+String((long) res.neuronsperline)+'\n';
  return res;
}

Shape_Rectangle_Result network::shape_rectangle(Command_Line_Parameters & clp) {
// gives neurons in the network positions such that they form a rectangle and
// are evenly distributed.
// horizlength and vertlength are the horizontal and vertical length of the rectangle in micrometers
  Shape_Rectangle_Result res;
  double horizlength = 10.0, vertlength = 10.0;
  int n;
  if ((n=clp.Specifies_Parameter("shape_horizontal"))>=0) horizlength = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("shape_vertical"))>=0) vertlength = atof(clp.ParValue(n));
  res.displaywidth = horizlength;
  int numneurons = PLLRoot<neuron>::length();
  if (numneurons<1) return res;
  netinfo += "shape rectangle with "+String(horizlength,"%.2f")+"horizontal and "+String(vertlength,"%.2f")+"vertical micrometer sides\n";
  res.rectangulararea = horizlength*vertlength;
  netinfo += "area "+String(res.rectangulararea,"%.3f")+" micrometers^2\n";
  double areaperneuron = res.rectangulararea/((double) numneurons);
  double squaresidelength = sqrt(areaperneuron);
  res.neuronsperline = (int) (horizlength/squaresidelength);
  if (res.neuronsperline<1) res.neuronsperline = 1;
  // NOTE: After the (int) approximation above, preceding calculations are
  // approximations! For instance, we do not use squaresidelength directly
  // for the neuron spacing below for that reason.
  int discardneurons = numneurons % res.neuronsperline;
  if (discardneurons>0) {
    for (int i = 0; i<discardneurons; i++) PLLRoot<neuron>::tail()->remove();
    numneurons = PLLRoot<neuron>::length();
    res.message = "Number of neurons reduced for even distribution.\nNew number of neurons = " + String((long) numneurons) + "\n";
  }
  netinfo += "shaped network consists of "+String((long) numneurons)+" neurons\n";
  double rectangleneuronxspacing = horizlength/((double) (res.neuronsperline-1)); // see reason in NOTE above
  double rectangleneuronyspacing = vertlength/((double) ((numneurons/res.neuronsperline)-1)); // see reason in NOTE above
  netinfo += "horizontal distance between neurons is "+String(rectangleneuronxspacing,"%.3f")+" micrometers\n";
  netinfo += "vertical distance between neurons is "+String(rectangleneuronyspacing,"%.3f")+" micrometers\n";
  int colnum = 0; double liney = 0.0;
  bool perturbpositions = false; double perturbfactorx = 0.0, perturbfactory =0.0;
  if ((paramspreadmin) && (paramspreadmax)) if (paramspreadmin->Pos().X()<paramspreadmax->Pos().X()) {
  	perturbpositions = true;
	perturbfactorx = (paramspreadmax->Pos().X()-paramspreadmin->Pos().X())/((double) X_misc.get_max());
	perturbfactory = (paramspreadmax->Pos().Y()-paramspreadmin->Pos().Y())/((double) X_misc.get_max());
     netinfo += "neuron positions perturbed horizontally between " + String(paramspreadmin->Pos().X(),"%.2f") + " and " + String(paramspreadmax->Pos().X(),"%.2f") + " micrometers\nneuron positions perturbed vertically between " + String(paramspreadmin->Pos().Y(),"%.2f") + " and " + String(paramspreadmax->Pos().Y(),"%.2f") + " micrometers\n";
  }
  double neuronx = 0.0, neurony = 0.0, perturbx, perturby;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    neurony = liney;
    neuronx = ((double) colnum)*rectangleneuronxspacing;
    if (perturbpositions) {
      perturbx = (perturbfactorx*((double) X_misc.get_rand())) + paramspreadmin->Pos().X();
      perturby = (perturbfactory*((double) X_misc.get_rand())) + paramspreadmin->Pos().Y();
	 e->set_position_in_Z_plane(neuronx+perturbx,neurony+perturby);
    } else e->set_position_in_Z_plane(neuronx,neurony);
    colnum++;
    if (colnum>=res.neuronsperline) {
      colnum=0;
      liney += rectangleneuronyspacing;
    }
  }
  center.set_X(neuronx/2.0);
  center.set_Y(neurony/2.0);
  // 3. insure that no more than one neuron occupies one position
  // Allocate to one network region
  neuron * rootneuron = PLLRoot<neuron>::head();
  region & reg = regions.append("net",*rootneuron);
  if (rootneuron) PLL_LOOP_FORWARD(neuron,rootneuron->Next(),1) reg.append(*e);
  res.message += "Number of neurons per line = "+String((long) res.neuronsperline)+'\n';
  return res;
}

std::vector<String> get_chem_factors(String factorsstr) {
    std::vector<String> chemfactors_vec;
    while (!factorsstr.empty()) {
      String factor(factorsstr.before(' '));
      if (factor.empty()) {
        factor = factorsstr;
        factorsstr = "";
      } else {
        factorsstr = factorsstr.after(' ');
      }
      chemfactors_vec.emplace_back(factor);
    }
    return chemfactors_vec;
}

#ifdef VECTOR3D
Shape_Box_Result network::shape_box(Command_Line_Parameters & clp) {
// gives neurons in the network positions such that they form a box and
// are evenly distributed.
// horizlength, vertlength and depth are the three dimensional sizes of the
// box in micrometers
  Shape_Box_Result res;
  double horizlength = 10.0, vertlength = 10.0, depth = 10.0;
  int n;
  if ((n=clp.Specifies_Parameter("shape_horizontal"))>=0) horizlength = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("shape_vertical"))>=0) vertlength = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("shape_depth"))>=0) depth = atof(clp.ParValue(n));
  res.displaywidth = horizlength;
  int numneurons = PLLRoot<neuron>::length();
  if (numneurons<1) return res;
  netinfo += "shape box with "+String(horizlength,"%.2f")+" horizontal, "+String(vertlength,"%.2f")+" vertical and "+String(depth,"%.2f")+" depth micrometer sides\n";
  res.boxvolume = horizlength*vertlength*depth;
  //[***REMOVE]res.rectangulararea = horizlength*vertlength;
  netinfo += "volume "+String(res.boxvolume,"%.3f")+" micrometers^3\n";
  double volumeperneuron = res.boxvolume/((double) numneurons);
  double cubesidelength = exp((1.0/3.0)*log(volumeperneuron)); // the cubed root
  res.neuronsperline = (int) (horizlength/cubesidelength);
  if (res.neuronsperline<1) res.neuronsperline = 1;
  res.neuronsdeep = (int) (depth/cubesidelength);
  if (res.neuronsdeep<1) res.neuronsdeep = 1;
  res.neuronsperxzplane = res.neuronsperline * res.neuronsdeep;
  // NOTE: After the (int) approximation above, preceding calculations are
  // approximations! For instance, we do not use cubesidelength directly
  // for the neuron spacing below for that reason.
  int discardneurons = numneurons % res.neuronsperxzplane;
  if (discardneurons>0) {
    for (int i = 0; i<discardneurons; i++) PLLRoot<neuron>::tail()->remove();
    numneurons = PLLRoot<neuron>::length();
    res.message = "Number of neurons reduced for even distribution.\nNew number of neurons = " + String((long) numneurons) + "\n";
  }
  netinfo += "shaped network consists of "+String((long) numneurons)+" neurons\n";
  double boxneuronxspacing = horizlength/((double) (res.neuronsperline-1)); // see reason in NOTE above
  double boxneuronzspacing = depth/((double) (res.neuronsdeep-1));
  double boxneuronyspacing = vertlength/((double) ((numneurons/res.neuronsperxzplane)-1)); // see reason in NOTE above
  netinfo += "horizontal distance between neurons is "+String(boxneuronxspacing,"%.3f")+" micrometers\n";
  netinfo += "vertical distance between neurons is "+String(boxneuronyspacing,"%.3f")+" micrometers\n";
  netinfo += "depth distance between neurons is "+String(boxneuronzspacing,"%.3f")+" micrometers\n";
  int colnum = 0, depthnum = 0; double liney = 0.0, linez = 0.0;
  bool perturbpositions = false;
  spatial perturbfactor;
  if ((paramspreadmin) && (paramspreadmax)) if (paramspreadmin->Pos().X()<paramspreadmax->Pos().X()) {
  	perturbpositions = true;
	perturbfactor.set_all((paramspreadmax->Pos().X()-paramspreadmin->Pos().X())/((double) X_misc.get_max()),(paramspreadmax->Pos().Y()-paramspreadmin->Pos().Y())/((double) X_misc.get_max()),(paramspreadmax->Pos().Z()-paramspreadmin->Pos().Z())/((double) X_misc.get_max()));
     netinfo += "neuron positions perturbed horizontally between " + String(paramspreadmin->Pos().X(),"%.2f") + " and " + String(paramspreadmax->Pos().X(),"%.2f") + " micrometers\nneuron positions perturbed vertically between " + String(paramspreadmin->Pos().Y(),"%.2f") + " and " + String(paramspreadmax->Pos().Y(),"%.2f") + " micrometers\nneuron positions perturbed in depth between " + String(paramspreadmin->Pos().Z(),"%.2f") + " and " + String(paramspreadmax->Pos().Z(),"%.2f") + " micrometers\n";
  }
  spatial npos, perturb;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    npos.set_all(((double) colnum)*boxneuronxspacing,liney,linez);
    if (perturbpositions) {
      perturb.set_all((perturbfactor.X()*((double) X_misc.get_rand())) + paramspreadmin->Pos().X(),(perturbfactor.Y()*((double) X_misc.get_rand())) + paramspreadmin->Pos().Y(),(perturbfactor.Z()*((double) X_misc.get_rand())) + paramspreadmin->Pos().Z());
      npos += perturb;
    }
    e->set_position(npos);
    colnum++;
    if (colnum>=res.neuronsperline) {
      colnum = 0;
      linez += boxneuronzspacing;
      depthnum++;
      if (depthnum>=res.neuronsdeep) {
	linez = 0.0;
	depthnum = 0;
	liney += boxneuronyspacing;
      }
    }
  }
  npos /= 2.0; center = npos;
  // 3. insure that no more than one neuron occupies one position
  // Allocate to one network region
  neuron * rootneuron = PLLRoot<neuron>::head();
  region & reg = regions.append("net",*rootneuron);
  if (rootneuron) PLL_LOOP_FORWARD(neuron,rootneuron->Next(),1) reg.append(*e);
  res.message += "Number of neurons per line = "+String((long) res.neuronsperline)+'\n';
  res.message += "Number of neurons deep     = "+String((long) res.neuronsdeep)+'\n';
  return res;
}

void disc_shaped_volume::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  String idstr("shape.");
  if (label) {
    if ((n=clp.Specifies_Parameter(*label+".centerX"))>=0) center.set_X(atof(clp.ParValue(n)));
    if ((n=clp.Specifies_Parameter(*label+".centerY"))>=0) center.set_Y(atof(clp.ParValue(n)));
    if ((n=clp.Specifies_Parameter(*label+".centerZ"))>=0) center.set_Z(atof(clp.ParValue(n)));
    idstr.prepend(*label+'.');
  }
  if ((n=clp.Specifies_Parameter(idstr+"radius"))>=0) radius = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(idstr+"thickness"))>=0) thickness = atof(clp.ParValue(n));
}

String disc_shaped_volume::report_parameters() {
  String res;
  if (label) res = *label + ' ';
  res += "disc shape: center = ("+String(center.X(),"%.3f,")+String(center.Y(),"%.3f,")+String(center.Z(),"%.3f) radius = ")+String(radius,"%.3f thickness = ")+String(thickness,"%.3f\n");
  return res;
}

spatial disc_shaped_volume::random_location() {
  // Returns a random location within the disc.
  /* spatial res(radius*random_double(),2.0*M_PI*random_double(),M_PI/2.0); [***NOTE] This creates a distribution clustered around the center!
  res.convert_from_spherical();
  double zperturbation = thickness/2.0;
  res.add_Z(random_double(-zperturbation,zperturbation));
  res += center; */
  double thhalf = thickness/2.0;
  while (1) {
    spatial res(X_misc.get_rand_range_real1(-radius,radius),X_misc.get_rand_range_real1(-radius,radius),X_misc.get_rand_range_real1(-thhalf,thhalf));
    spatial testin(res); testin.set_Z(0.0); // flatten testin into a plane
    testin.convert_to_spherical();
    if (testin.X()<=radius) {
      res += center;
      return res;
    }
  }
  return spatial();
}

void box_shaped_volume::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  String idstr("shape.");
  if (label) {
    if ((n=clp.Specifies_Parameter(*label+".centerX"))>=0) center.set_X(atof(clp.ParValue(n)));
    if ((n=clp.Specifies_Parameter(*label+".centerY"))>=0) center.set_Y(atof(clp.ParValue(n)));
    if ((n=clp.Specifies_Parameter(*label+".centerZ"))>=0) center.set_Z(atof(clp.ParValue(n)));
    idstr.prepend(*label+'.');
  }
  if ((n=clp.Specifies_Parameter(idstr+"height"))>=0) height = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(idstr+"width"))>=0) width = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(idstr+"depth"))>=0) depth = atof(clp.ParValue(n));
}

String box_shaped_volume::report_parameters() {
  String res;
  if (label) res = *label + ' ';
  res += "box shape: center = ("+String(center.X(),"%.3f,")+String(center.Y(),"%.3f,")+String(center.Z(),"%.3f) width = ")+String(width,"%.3f height = ")+String(height,"%.3f depth = ")+String(depth,"%.3f\n");
  return res;
}

spatial box_shaped_volume::random_location() {
  // Returns a random location within the box.
  spatial res(width*X_misc.get_rand_real1(),height*X_misc.get_rand_real1(),depth*X_misc.get_rand_real1());
  spatial boxmiddlevec(width/2.0,height/2.0,depth/2.0);
  res -= boxmiddlevec;
  res += center;
  return res;
}

void sphere_shaped_volume::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  String idstr("shape.");
  if (label) {
    if ((n=clp.Specifies_Parameter(*label+".centerX"))>=0) center.set_X(atof(clp.ParValue(n)));
    if ((n=clp.Specifies_Parameter(*label+".centerY"))>=0) center.set_Y(atof(clp.ParValue(n)));
    if ((n=clp.Specifies_Parameter(*label+".centerZ"))>=0) center.set_Z(atof(clp.ParValue(n)));
    idstr.prepend(*label+'.');
  }
  if ((n=clp.Specifies_Parameter(idstr+"radius"))>=0) radius = atof(clp.ParValue(n));
}

String sphere_shaped_volume::report_parameters() {
  String res;
  if (label) res = *label + ' ';
  res += "sphere shape: center = ("+String(center.X(),"%.3f,")+String(center.Y(),"%.3f,")+String(center.Z(),"%.3f) radius = ")+String(radius,"%.3f\n");
  return res;
}

spatial sphere_shaped_volume::random_location() {
  // Returns a random location within the sphere.
  //spatial res(radius*X_misc.get_rand_real1(),2.0*M_PI*X_misc.get_rand_real1(),M_PI*X_misc.get_rand_real1()); // *** This way favors locations in the center.
  spatial res;
  double diameter = 2.0*radius;
  double radiusSQ = radius*radius;
  do {
    res.set_all((diameter*X_misc.get_rand_real1())-radius,(diameter*X_misc.get_rand_real1())-radius,(diameter*X_misc.get_rand_real1())-radius);
  } while (res.len2()>radiusSQ);
  res += center;
  return res;
}

region_parameters::region_parameters(const char * l, Command_Line_Parameters & clp): label(l), sv(NULL), memberneurons(NULL), nummemberneurons(0), generalandspecificneurons(0), minneuronseparation(10.0) {
  int n;
  if ((n=clp.Specifies_Parameter(label+".shape"))>=0) {
    String shape(clp.ParValue(n));
    if (shape=="disc") {
      sv = new disc_shaped_volume(&label);
      sv->parse_CLP(clp);
    } else if (shape=="box") {
      sv = new box_shaped_volume(&label);
      sv->parse_CLP(clp);
    } else if (shape=="sphere") {
      sv = new sphere_shaped_volume(&label);
      sv->parse_CLP(clp);
    } else error("Unknown region shape ("+shape+").\n");
  }
  if ((n=clp.Specifies_Parameter(label+".neurons"))>=0) nummemberneurons = atoi(clp.ParValue(n));
  generalandspecificneurons = nummemberneurons;
  for (neuron_type nt = neuron_type(0); nt < UNTYPED_NEURON; nt = neuron_type(nt+1)) { // Read the number of neurons of a specific type to be added to the region.
    if ((n=clp.Specifies_Parameter((label+'.')+neuron_short_name[nt]))>=0) {
      specificneurons[nt] = atoi(clp.ParValue(n));
      //cout << label << '.' << neuron_short_name[nt] << '=' << specificneurons[nt] << '\n'; cout.flush();
    } else specificneurons[nt] = 0;
    generalandspecificneurons += specificneurons[nt];
  }
  if ((n=clp.Specifies_Parameter(label+".minneuronseparation"))>=0) minneuronseparation = atof(clp.ParValue(n));
  if (!sv) {
    warning("Warning: Region "+label+" has no specific shape, defaulting to disc.\n");
    sv = new disc_shaped_volume(&label);
    sv->parse_CLP(clp);
  }
}

String region_parameters::report_parameters() {
  String res((' '+String((long) nummemberneurons))+" general pool neurons");
  for (neuron_type nt = neuron_type(0); nt < UNTYPED_NEURON; nt = neuron_type(nt+1))
    if (specificneurons[nt]>0) res += (((" + "+String((long) specificneurons[nt]))+' ')+neuron_short_name[nt]) + " cells";
  return res;
}

#ifdef TESTING_SIMPLE_ATTRACTION
int attracts = 0;
#endif

neuron * region_parameters::add_neurons(PLLRoot<neuron> * all, neuron * n, neuron * psmin, neuron * psmax) {
  // The requisite number of neurons are recruited from the available
  // list in n and are given positions in the region.
  // Also, specific neuron types allocated to the region are created,
  // initialized and added to the master list of neurons. (See TL#200807040256.1.)
  if ((generalandspecificneurons<1) || (!all)) {
    warning("Warning: Region "+label+" is not associated with any neurons.\n");
    return n;
  }
  if (memberneurons) {
    warning("Warning: Region "+label+" attempted to add_neurons() more than once.\n");
    return n;
  }
  if (nummemberneurons) if ((!n) || (n->length()<nummemberneurons)) {
    int n_num = 0;
    if (n) n_num = n->length();
    generalandspecificneurons -= (nummemberneurons-n_num); // remove what cannot be allocated
    nummemberneurons = n_num;
    warning("Warning: Insufficient general pool neurons available for region "+label+".\n");
  }
  memberneurons = new neuronptr[generalandspecificneurons];
  int i = 0, attempts;
  PLL_LOOP_FORWARD(neuron,n,i<nummemberneurons) memberneurons[i++] = e;
  neuron * generalnext;
  if (nummemberneurons<1) generalnext = n;
  else generalnext = memberneurons[nummemberneurons-1]->Next(); // store this, just in case
  // Here add specific neurons and set their initial parameters as was done for the general pool neurons
  for (neuron_type nt = neuron_type(0); nt < UNTYPED_NEURON; nt = neuron_type(nt+1)) {
    for (int j = 0; j<specificneurons[nt]; j++) {
      // create and add
      switch (nt) {
      case PRINCIPAL_NEURON: memberneurons[i] = new principal(); break;;
      case INTERNEURON: memberneurons[i] = new interneuron(); break;;
      case MULTIPOLAR_NONPYRAMIDAL: memberneurons[i] = new multipolar_nonpyramidal(); break;;
      case BIPOLAR: memberneurons[i] = new bipolar(); break;;
      case PYRAMIDAL: memberneurons[i] = new pyramidal(); break;;
      default: memberneurons[i] = new principal(); break;; //neuron(); break;;
      }
      memberneurons[i]->parse_CLP(*main_clp); // This calls parameters parsing on neuron::parse_CLP() (a virtual function).
      // initialize
      if ((psmin) && (psmax)) if (psmin->Radius()<psmax->Radius()) {
	memberneurons[i]->set_radius(X_misc.get_rand_range_real1(psmin->Radius(),psmax->Radius()));
      }
      // add to the master list of neurons maintained by the network object (prepend to avoid interfering with parameter n)
      all->link_after(memberneurons[i]);
      i++;
    }
  }
#ifdef TESTING_SIMPLE_ATTRACTION
  if (eq) {
    network * net_ptr = eq->Net();
    int clp_n;
    if ((clp_n=main_clp->Specifies_Parameter(label+".attractors"))>=0) {
      // separate list of attractor factors
      std::vector<String> factors = get_chem_factors(main_clp->URI_unescape_ParValue(clp_n));
      for (auto & factor : factors) {
        // the neuron index by attractor reference in the network list
        auto chemfactor = net_ptr->get_or_add_chemlabel(factor);
        for (i =0; i<generalandspecificneurons; i++) {
          net_ptr->chemdata.attractor_somata[chemfactor].emplace(memberneurons[i]);
          memberneurons[i]->chemdata.attractor.emplace(chemfactor);
          net_ptr->target_neurons.emplace(memberneurons[i]);
        }
      }
      net_ptr->chemdata.has_specified_factors = true;
    }
    if ((clp_n=main_clp->Specifies_Parameter(label+".attractedto"))>=0) {
      // separate list of attractor factors
      std::vector<String> factors = get_chem_factors(main_clp->URI_unescape_ParValue(clp_n));
      for (auto & factor : factors) {
        // store them for each neuron in the region
        auto chemfactor = net_ptr->get_or_add_chemlabel(factor);
        for (i =0; i<generalandspecificneurons; i++) {
          memberneurons[i]->chemdata.attractedto.emplace(chemfactor);
        }
      }
      net_ptr->chemdata.has_specified_factors = true;
    }

    if (!net_ptr->chemdata.has_specified_factors) {
      // This test was meant for 2 layers that are attracted to each other.
      for (i = 0; i<generalandspecificneurons; i++) {
        memberneurons[i]->attracts = attracts;
        memberneurons[i]->attractedto = 1 - attracts;
      }
    }
  }
#endif
  // Set neuron positions.
  double SQminsep = minneuronseparation*minneuronseparation, sep, cellextents;
  for (i = 0; i<generalandspecificneurons; i++) {
    for (attempts = 1000; attempts>0; attempts--) {
      spatial npos(sv->random_location());
      // Check distance to other neurons.
      bool validlocation = true;
      PLL_LOOP_FORWARD(neuron,all->head(),e!=generalnext) if (e!=memberneurons[i]) { // Check (e!=generalnext) in order not to test with unplaced general pool neurons.
	if ((sep=distance2(e->Pos(),npos))<SQminsep) {
	  validlocation = false;
	  break;
	} else {
	  // Additional test in case cells need more space than
	  // minneuronseparation.
	  cellextents = e->Radius() + memberneurons[i]->Radius();
	  cellextents *= cellextents;
	  if (sep<cellextents) {
	    validlocation = false;
	    break;
	  }
	}
      }
      if (validlocation) {
	memberneurons[i]->set_position(npos);
	break;
      }
    }
    if (attempts<=0) error("Unable to find room for another neuron in region "+label+" after 1000 attempts\n");
  }
  return generalnext;
}

double region_parameters::average_distance_to_nearest_neighbor() {
  if (!memberneurons) {
    warning("Warning: Average distance to nearest neighbor reqeusted for region with no member neurons.\n");
    return -1.0;
  }
  double avdistnearest = 0.0;
  for (int i = 0; i<nummemberneurons; i++) {
    double nearest = MAXDOUBLE, newnearest;
    for (int j = 0; j<generalandspecificneurons; j++) if (i!=j) if ((newnearest=distance2(memberneurons[i]->Pos(),memberneurons[j]->Pos()))<nearest) nearest = newnearest;
    avdistnearest += sqrt(nearest);
  }
  return avdistnearest /= (double) generalandspecificneurons;
}

void region_parameters_root::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("regions"))>=0) {
    String regionstr(clp.URI_unescape_ParValue(n));
    while (!regionstr.empty()) {
      String regionname(regionstr.before(' '));
      if (regionname.empty()) {
        regionname = regionstr;
        regionstr = "";
      } else {
        regionstr = regionstr.after(' ');
      }
      link_before(new region_parameters(regionname,clp));
    }
  }
  numregions = length();
  if (numregions<1) error("Invalid number of regions ("+String((long) numregions)+").\n");
}

String region_parameters_root::report_parameters() {
  String res;
  if (numregions>0) {
    res = "Regions:\n";
    PLL_LOOP_FORWARD(region_parameters,head(),1) res += "  " + e->Label() + " (" + e->shape().str() + "):" + e->report_parameters() + '\n';
  }
  return res;
}

Shape_Regions_Result network::shape_regions(Command_Line_Parameters & clp,neuron * psmin,neuron * psmax) {
// Distributes neurons randomly in network positions according to
// specified regions. This is useful for the generation of cortical layers.
// Neurons are evenly distributed.
// Sizes are indicated in micrometers.
  Shape_Regions_Result res;
  // [***INCOMPLETE] This is a very rudimentary implementation with the sole
  // purpose of allowing me to explore the growth of networks in layers.
  // *** Possibly move the next three lines to nibr.cc.
  region_parameters_root rp;
  rp.parse_CLP(clp);
  res.message = rp.report_parameters();
  if ((numneurons_is_default) && (rp.total_general_pool_N()>0)) warning("Warning: Regions expect general pool neurons, but no neurons were assigned to the general pool throught the 'neurons' command with approximate proportions or through exact population numbers. Attempting to allocate general pool neurons from a general pool of default size.\n");
  if (rp.total_N()<1) {
    warning("Warning: Regions were defined, but no general pool or neuron type specific populations of neurons specified. A default of one general pool neuron will be specified for each region.");
    PLL_LOOP_FORWARD(region_parameters,rp.head(),1) e->set_general_pool(1);
    if (numneurons_is_default) warning("Warning: Allocating general pool neurons from a general pool of default size.\n");
  }
  neuron * nextavailable = PLLRoot<neuron>::head();
  double avnearestdist = 0.0;
  PLL_LOOP_FORWARD(region_parameters,rp.head(),1) {
    netinfo += e->shape().report_parameters();
    if (e->shape().displaywidth()>res.displaywidth) res.displaywidth = e->shape().displaywidth();
    // Feed the neurons into the region
    nextavailable = e->add_neurons(this,nextavailable,psmin,psmax);
#ifdef TESTING_SIMPLE_ATTRACTION
    attracts++; // For this test, the attracts label is simply increased for each region.
#endif
    if (e->N()>0) {
      avnearestdist += e->average_distance_to_nearest_neighbor();
      // Allocate to one network region
      neuronptr * nptrlist = e->neurons();
      region & reg = regions.append(e->Label(),*(nptrlist[0]));
      for (int i = 1; i<e->N(); i++) reg.append(*(nptrlist[i]));
    }
  }
  if (nextavailable) { // There were more "general pool" neurons than the sum of all "general pool member neurons" in the specified regions
    if (!numneurons_is_default) warning("Warning: More general pool neurons were provided to the network than the sum of general pool member neurons specified in all regions. The remaining general pool neurons will be discarded.\n");
    neuron * next_n_remove = nextavailable;
    for (neuron * n_remove = nextavailable; (n_remove); n_remove = next_n_remove) {
      next_n_remove = n_remove->Next();
      n_remove->remove();
    }
  }
  avnearestdist /= (double) rp.length();
  netinfo += "average distance to the nearest neighbor neuron within a region is "+String(avnearestdist,"%.3f\n");
  return res;
}
#endif // VECTOR3D

Shape_Circle_Result network::shape_circle(Command_Line_Parameters & clp) {
// gives neurons in the network positions such that they form a circle and
// are evenly distributed.
// the radius of the circle is given in micrometers
  Shape_Circle_Result res;
  double radius = 10.0, minneuronseparation = 0.0; bool randompolar = true;
  int n;
  if ((n=clp.Specifies_Parameter("shape_horizontal"))>=0) radius = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("shape_radius"))>=0) radius = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("shape_randompolar"))>=0) randompolar = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("minneuronseparation"))>=0) minneuronseparation = atof(clp.ParValue(n));
  res.displaywidth = 2.0*radius;
  int numneurons = PLLRoot<neuron>::length();
  if (numneurons<1) return res;
  netinfo += "shape cicle with "+String(radius,"%.2f")+" micrometer radius\n";
  res.circlearea = M_PI*radius*radius;
  netinfo += "area "+String(res.circlearea,"%.3f")+" micrometers^2\n";
  // the circle is placed within a square with sides equal to the diameter
  // and neurons beyond the circle boundary are discarded
  double areaperneuron = (4.0*radius*radius)/((double) numneurons);
  double squaresidelength = sqrt(areaperneuron);
  res.maxneuronsperline = (int) ((2.0*radius)/squaresidelength);
  if (res.maxneuronsperline<1) res.maxneuronsperline = 1;
  // NOTE: After the (int) approximation above, preceding calculations are
  // approximations! For instance, we do not use squaresidelength directly
  // for the neuron spacing below for that reason.
  int removenum = 0;
  if (!randompolar) {
    int discardneurons = numneurons % res.maxneuronsperline;
    if (discardneurons>0) {
      for (int i = 0; i<discardneurons; i++) PLLRoot<neuron>::tail()->remove();
      numneurons = PLLRoot<neuron>::length();
      res.message = "Number of neurons reduced for even distribution.\nNew number of neurons = " + String((long) numneurons) + "\n";
    }
    double rectangleneuronxspacing = (2.0*radius)/((double) (res.maxneuronsperline-1)); // see reason in NOTE above
    double rectangleneuronyspacing = (2.0*radius)/((double) ((numneurons/res.maxneuronsperline)-1)); // see reason in NOTE above
    int colnum = 0; double liney = 0.0;
    double neuronx = 0.0, neurony = 0.0;
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
      neurony = liney;
      neuronx = ((double) colnum)*rectangleneuronxspacing;
      e->set_position_in_Z_plane(neuronx,neurony);
      colnum++;
      if (colnum>=res.maxneuronsperline) {
	colnum=0;
	liney += rectangleneuronyspacing;
      }
    }
    // remove neurons beyond the circle boundary
    double sqrradius = radius*radius;
    center.set_X(neuronx/2.0);
    center.set_Y(neurony/2.0);
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) if (distance2(center,e->Pos())>sqrradius) removenum++;
    if (numneurons>removenum) {
      bool removeprev = false;
      PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
	if (removeprev) e->Prev()->remove();
	if (distance2(center,e->Pos())>sqrradius) removeprev = true;
	else removeprev = false;
      }
      if (removeprev) PLLRoot<neuron>::tail()->remove();
      numneurons = PLLRoot<neuron>::length();
      res.message = "Number of neurons reduced for circle shape.\nNew number of neurons = " + String((long) numneurons) + "\n";
      netinfo += "horizontal distance between neurons is "+String(rectangleneuronxspacing,"%.3f")+" micrometers\n";
      netinfo += "vertical distance between neurons is "+String(rectangleneuronyspacing,"%.3f")+" micrometers\n";
      // 3. insure that no more than one neuron occupies one position
      res.message += "Maximum number of neurons per line = "+String((long) res.maxneuronsperline)+'\n';
    } else res.message = "Switching to random polar placement to avoid a reduction to zero neurons\n";
  }
  if ((randompolar) || (numneurons<=removenum)) {
    // Setting neuron positions randomly in polar coordinates
    //double sumx = 0.0, sumy = 0.0;
    int attempts;
    double SQminsep = minneuronseparation*minneuronseparation, sep, cellextents, x, y, sqrradius = radius*radius;
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
      // Set neuron positions (as in regions shapes).
      for (attempts = 1000; attempts>0; attempts--) {
	//double r = radius*X_misc.get_rand_real1(); // *** can't do this, increases density toward center
	//double theta = 2.0*M_PI*X_misc.get_rand_real1(); // *** can't do this, increases density toward center
	//e->set_position_in_Z_plane(r,theta); // *** can't do this, increases density toward center
	//e->Pos().convert_from_spherical(); // *** can't do this, increases density toward center
	do { // allow only random locations within the circle
	  x = radius*X_misc.get_rand_range_real1(-1.0,1.0);
	  y = radius*X_misc.get_rand_range_real1(-1.0,1.0);
	} while (((x*x) + (y*y)) > sqrradius);
	e->set_position_in_Z_plane(x,y);
	bool validlocation = true;
	PLL_LOOP_FORWARD_NESTED(neuron,PLLRoot<neuron>::head(),(ne!=e),ne) {
	  if ((sep=distance2(ne->Pos(),e->Pos()))<SQminsep) {
	    validlocation = false;
	    break;
	  } else {
	    // Additional test in case cells need more space than
	    // minneuronseparation.
	    cellextents = ne->Radius() + e->Radius();
	    cellextents *= cellextents;
	    if (sep<cellextents) {
	      validlocation = false;
	      break;
	    }
	  }
	}
	if (validlocation) break;
      }
      if (attempts<=0) error("Unable to find room for another neuron in circle after 1000 attempts\n");
      //sumx += e->Pos().X(); sumy += e->Pos().Y();
      // Check distance to other neurons.
    }
    numneurons = PLLRoot<neuron>::length();
    center.set_X(0.0); // sumx/(double) numneurons
    center.set_Y(0.0); // sumy/(double) numneurons
  }
  // Allocate to one network region
  neuron * rootneuron = PLLRoot<neuron>::head();
  region & reg = regions.append("net",*rootneuron);
  if (rootneuron) PLL_LOOP_FORWARD(neuron,rootneuron->Next(),1) reg.append(*e);
  netinfo += "shaped network consists of "+String((long) numneurons)+" neurons\n";
  return res;
}

Shape_Result network::shape_network(Command_Line_Parameters & clp, neuron * psmin, neuron * psmax) {
  String shapename("rectangle");
  int n;
  if ((n=clp.Specifies_Parameter("shape"))>=0) shapename = downcase(clp.ParValue(n));
  if (shapename=="rectangle")  return shape_rectangle(clp);
#ifdef VECTOR3D
  else if (shapename=="box") return shape_box(clp);
  else if (shapename=="regions") return shape_regions(clp,psmin,psmax);
#endif
  else if (shapename=="circle") return shape_circle(clp);
  else if (shapename=="hexagon") return shape_hexagon(clp);
  error("Error: Unrecognized shape.\n");
  return Shape_Result();
}

neuron * network::find_neuron_by_idx(long idx) {
  PLL_LOOP_FORWARD(neuron, PLLRoot<neuron>::head(), 1) {
    if (e->numerical_ID()==idx) return e;
  }
  return NULL;
}

neuron * network::nearest(spatial & p) {
  neuron * res = PLLRoot<neuron>::head();
  if (!res) return NULL;
  double d, mindist = distance2(p,PLLRoot<neuron>::head()->Pos());
  PLL_LOOP_FORWARD(neuron,res->Next(),1) if ((d=distance2(p,e->Pos()))<mindist) {
    res = e;
    mindist = d;
  }
  return res;
}

neuronset * network::inrange(neuron * n, double radius, neuron_type ntype) {
  // gets a set of neurons that are in range of a radius
  // this may loop around edges if edges==false
  // ntype allows specification of a type of neuron to target
  spatial p(n->Pos());
  // [***INCOMPLETE] Perhaps there should be a 3D equivalent of the torus
  // connectivity that avoids edge effects.
#ifdef VECTOR2D
  spatial p2, p3;
  if (!edges) {
    if (p.Y()<center.Y()) p2.set_Y(p.Y()+(2.0*center.Y()));
    else p2.set_Y(p.Y()-(2.0*center.Y()));
    p2.set_X((2.0*center.X())-p.X());
    if (p.X()<center.X()) p3.set_X(p.X()+(2.0*center.X()));
    else p3.set_X(p.X()-(2.0*center.X()));
    p3.set_Y((2.0*center.Y())-p.Y());
  }
#endif
  radius *= radius;
  neuronset * nset = new neuronset();
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) if ((e!=n) && ((ntype==UNTYPED_NEURON) || (ntype==e->TypeID()))) {
    if (distance2(e->Pos(),p)<=radius) nset->link_before(new neuronsetel(*e));
#ifdef VECTOR2D
    else if (!edges) { // check if is in range when connected as a torus
      if (distance2(e->Pos(),p2)<=radius) nset->link_before(new neuronsetel(*e));
      else {
	if (distance2(e->Pos(),p3)<=radius) nset->link_before(new neuronsetel(*e));
      }
    }
#endif
  }
  return nset;
}

void network::uniform_random_connectivity(double range, int minconn, int maxconn, neuron_type sourcetype, neuron_type targettype, bool clearexisting) {
  // selects a set of neurons within the range around each neuron and creates a
  // random number (between minconn and maxconn) of connections between the
  // neurons and a uniform random selection of neurons from the set.
  // sourcetype and targettype allow selection of specific source and target
  // neuron types for which to create connectivity
  int numsynapses; int progresscounter = 0, progressfactor = 0;
  netinfo += "uniform random connectivity, random seed = " + String((long) random_seed) + "\nconnection range " + String(range,"%.1f") + ", with " + String((long) minconn) + " to " + String((long) maxconn) + " connections per neuron\n";
  if (outattr_show_progress) {
    progresscounter = progressfactor = (int) (((double) PLLRoot<neuron>::length())/25.0);
    if (progressfactor!=0) progress("progress: ");
  }
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) if ((sourcetype==UNTYPED_NEURON) || (sourcetype==e->TypeID())) {
    if (progressfactor!=0) {
      if ((progresscounter--)<=0) {
	progress("+");
	progresscounter = progressfactor;
      }
    }
    if (maxconn==minconn) numsynapses = maxconn;
    else numsynapses = minconn + (int) (((double) X_misc.get_rand())/((double) X_misc.get_max())*((double) (maxconn-minconn)));
    neuronset * nset = inrange(e,range,targettype);
    if (!nset) error("nibr Error: Empty neuron set in network::uniform_random_connectivity()\n");
    int idxlen = nset->length();
    // clear pre-existing connections, then connect to all in set, array to connections (rapid access)
    if (clearexisting) e->clear_output_connections();
    connectionptr * cset = new connectionptr[idxlen];
    int i = 0;;
    PLL_LOOP_FORWARD_NESTED(neuronsetel,nset->head(),1,t) {
      if (i>=idxlen) error("nibr Error: Attempting to make too many connecitons in network::uniform_random_connectivity()\n");
      cset[i] = e->connect_to(t->Neuron());
      i++;
    }
    delete nset;
    int nidx, syntype;
    // set up uniform random connectivity
    for (i = 0; i < numsynapses; i++) {
      nidx = X_misc.get_rand() % idxlen; // *** is this a valid pseudo-random method?
      // this is done here so that synapse type selection can also be
      // randomized - and the actual building of the synapse and of
      // possible arborization and morphology is done within the synapse
      // object
      // [***INCOMPLETE] syntype = possible_syntypes(cset[nidx]->PreSynaptic(),cset[nidx]->PostSynaptic());
      switch (syntype) {
      default: new candidate_synapse(*(cset[nidx]),1.0);
      }
    }
    //    cout << "#" << String((long) e->number_of_synapses());
    // remove connections without synapses
    for (i = 0; i < idxlen; i++) if (!cset[i]->Synapses()->head()) cset[i]->remove();
    delete[] cset;
  }
  if (progressfactor!=0) progress('\n');
  netinfo += "actual number of connections in network is " + String((long) number_of_connections()) + "\nactual number of synapses in network is " + String((long) number_of_synapses()) +'\n';
}

void network::develop_connection_structure(Connection_Statistics_Root & cstats, dendritic_growth_model & dgm, dendritic_growth_model & agm) {
  // <A NAME="network::develop_connection_structure"> </A>
  // This network interface supports methods for the generation of
  // presynaptic and postsynaptic structure that involve simulated
  // development (growth).
  netinfo += "developmental connection structure, random seed = " + String((long) random_seed) + "\nConnection Statistics:\n";
  if (!cstats.head()) error("nibr Error: Missing Connection_Statistics in network::develop_connection_structure()\n");
  netinfo += cstats.head()->all_str();
  // (this is now a class variable) double t = 0.0; // time in seconds
  dgm.initialize(this,&cstats,t); // also sets t = _dt if FIXED_STEP_SIMULATION_START_AT_DT
  agm.initialize(this,&cstats,t);
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
  if ((!initiallengthatactualfirstbranch) && (t<dgm.get_fixed_time_step_size())) t = dgm.get_fixed_time_step_size(); // (See TL#200603160459.1.)
#endif
  if (dgm.get_fixed_time_step_size()!=agm.get_fixed_time_step_size()) warning("Warning: fixed step sizes for dendritic and axonal growth are not equal in network;:develop_connection_structure()\n");
  initialize_spatial_segment_subset(cstats,dgm,agm);
  // Note that sampled_output calls are made here, since sampling is done
  // at the same time points for the entire network structure.
  if (outattr_show_progress) progress("Growing dendrites and axons in time steps of "+String(dgm.get_fixed_time_step_size(),"%.1f")+" seconds\nup to t = "+String(max_growth_time,"%.1f")+" seconds:\n0%                                   50%                                 100%\n");
  double dt_mark = max_growth_time/78.0; // 50.0;
  double t_nextmark = dt_mark;
  while (netmorph_active) {
    //cout << "t=" << t << '\n';
    if (sampled_output) sampled_output->Sample(t);
    if (progress_percentage != nullptr) {
      //double pct = 100.0*(t/max_growth_time);
      double pct = 99.0*(t/max_growth_time); // Let's keep 1% for the embedded model builder.
      (*progress_percentage) = int(pct); // update external progress tracker (e.g. used by NES)
    }
    if (t>=t_nextmark) {
      if (outattr_show_progress && (!embedlog)) progress('=');
      t_nextmark += dt_mark;
    }
    dgm.grow(this,&cstats,t);
    agm.grow(this,&cstats,t);
    // [***INCOMPLETE]
    // *** For connectivity, apply the connection equations to axons, soma
    // and dendrites within a segment, as well as within segments that are
    // immediate neighbors (e.g. above is identified in the hierarchy as
    // any subdivision with a bottom value equal to the testing segment top).
    // Connection equations may be applied whenever a segment is added to
    // the partition hierarhcy or at the very end of the growth simulation,
    // or after each step. I can make that a command line option.
    if (synapses_during_development) sfm->establish_connections();
    if (dgm.IsComplete(t)) break;
    if (completion.is_complete()) break;
  }
  if (outattr_show_progress) progress("\ndone.\nPerforming growth post-op tasks.\n");
  dgm.postop(this,&cstats,t);
  if (outattr_show_progress) progress(" Dendrites post-op complete.\n");
  agm.postop(this,&cstats,t);
  if (outattr_show_progress) progress(" Axons post-op complete.\n");
#ifdef TEST_FOR_NAN
  progress(" TEST_FOR_NAN: Not performing synapses post-op.");
#else
  sfm->establish_connections();
  remove_abstract_connections_without_synapses();
  if (outattr_show_progress) progress(" Synapses post-op complete.\n");
#endif
  if (sampled_output) sampled_output->Sample(t);
  if (ndm) ndm->diameters();
  if (outattr_show_progress) progress(" Post-op complete.\n");
}

void network::spatial_extents(spatial & minpoint, spatial & maxpoint) {
  minpoint = PLLRoot<neuron>::head()->Pos(); minpoint.add_to_all(-PLLRoot<neuron>::head()->Radius());
  maxpoint = PLLRoot<neuron>::head()->Pos(); maxpoint.add_to_all(PLLRoot<neuron>::head()->Radius());
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head()->Next(),1) {
    spatial minp(e->Pos()), maxp(e->Pos());
    minp.add_to_all(-e->Radius()); maxp.add_to_all(e->Radius());
    for (int i=0; i<PDIMSIZE; i++) {
      if (minp[i]<minpoint[i]) minpoint[i] = minp[i];
      if (maxp[i]>maxpoint[i]) maxpoint[i] = maxp[i];
    }
  }
}

void network::initialize_spatial_segment_subset(Connection_Statistics_Root & cstats, dendritic_growth_model & dgm, dendritic_growth_model & agm) {
  delete sss;
  // [***NOTE] There may be a way to grow the dimensions of the outmost spatial
  // segment when network components begin to extend beyond it, which would add
  // the ability to subdivide the space. Currently, anything that extends
  // beyond the formal boundaries of the outmost segment is simply a member of
  // that segment. A dynamic resizing of the outmost segment may involve
  // propagating constraint changes to some of the subdividing segments.
  fibre_structure * fs = PLLRoot<neuron>::head()->OutputStructure()->head();
  if (!fs) error("nibr error: missing axonal arbor in network::initialize_spatial_segment_subset()\n");
  if (!fs->ElongationModel()) error("nibr error: missing elongation model in network::initialize_spatial_segment_subset()\n");
  // [***BEWARE] Do not allow asynchronous event handling to occur during this
  // glimpse ahead!
  // --- start glimpse
  double t_cache = t;
  t = max_growth_time;
  double mae = fs->ElongationModel()->elongation(fs) / ((double) fs->count_terminal_segments());
  fs->ElongationModel()->reset();
  t = t_cache;
  // --- end glimpse
  mae *= 2.0;
  if (specified_max_spatial_segment_coordinate>0.0) mae = specified_max_spatial_segment_coordinate;
  double proximity_threshold = 0.0;
  for (int i=0; i<=UNTYPED_NEURON; i++) for (int j=0; j<=UNTYPED_NEURON; j++) if (maxfibredistance[i][j]>proximity_threshold) proximity_threshold = maxfibredistance[i][j];
  double subminspan = 2.0*proximity_threshold;
  if (minimum_span_spatialsegmentsubset>subminspan) subminspan = minimum_span_spatialsegmentsubset;
  report("Spatial segment subset minimum span distance: "+String(subminspan,"%.3f")+" [micron]\n");
  spatial minpoint, maxpoint;
  spatial_extents(minpoint,maxpoint);
  minpoint.add_to_all(-mae);
  maxpoint.add_to_all(mae);
  sss = new Spatial_Segment_Subset(NULL,0,spatialsegmentsubsetsizelimit,subminspan,minpoint,maxpoint);
  sfm->set_proximity_threshold(proximity_threshold);
}

int network::number_of_connections() {
  int sum = 0;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) sum += e->OutputConnections()->length();
  return sum;
}

int network::number_of_synapses() {
  int sum = 0;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) sum += e->number_of_synapses();
  return sum;
}

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
void network::Collect_Data(network_statistics_base & nsb) {
  nsb.Set_Tail_Age(t);
  if (nsb.Collecting_Statistics()) PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->Collect_Data(nsb);
}
#endif

double network::mean_number_of_output_terminal_segments(network_statistics_data * nsd) {
  // nt is the mean number of terminal segments per axonal tree
  // (multipolar neurons have multiple axonal trees)
  // std can return the standard deviation from that mean
  // arborsamples can return the number of arbors over which the statistic
  // was collected (so that a measure of statistical significance is possible)
  double mean_nt, d_nt, sumsq_nt = 0.0, numtrees = 0.0, nt;
  int sum_nt = 0;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    numtrees += (double) e->OutputStructure()->length();
    sum_nt += e->total_output_terminal_segments();
  }
  mean_nt = ((double) sum_nt) / numtrees;
  if (nsd) {
    nsd->n = numtrees; // the number of arbor samples
    nsd->mean = mean_nt;
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1)
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,pss) {
      nt = ((double) (pss->count_terminal_segments()));
      if (nt<nsd->min) nsd->min = nt;
      if (nt>nsd->max) nsd->max = nt;
      d_nt = nt - mean_nt;
      sumsq_nt += (d_nt)*(d_nt);
    }
    if (numtrees==1.0) nsd->std = -1.0;
    else {
      sumsq_nt /= (numtrees - 1.0);
      nsd->std = sqrt(sumsq_nt); // standard deviation sqrt((sumsq(a-mean(a)))/(n-1))
    }
  }
  return mean_nt;
}

double network::mean_presynaptic_structure_length(network_statistics_data & termnsd, network_statistics_data & nsd, network_statistics_data * internsd) {
  // This returns the mean length of axonal arbor in the network and the
  // mean length of axonal terminal segments in the network.
  double arborlength, totterminalseglength, lendiff;
  int numtrees = 0;
  nsd.mean = 0.0;
  termnsd.n = 0.0;
  termnsd.mean = 0.0;
  if (internsd) {
    internsd->mean = 0.0;
    internsd->n = 0.0;
  }
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1)
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,pss) {
    arborlength = pss->total_arbor_length(&totterminalseglength,internsd); // also cached in pss
    if (arborlength<nsd.min) nsd.min = arborlength;
    if (arborlength>nsd.max) nsd.max = arborlength;
    nsd.mean += arborlength;
    termnsd.mean += totterminalseglength;
    termnsd.n += ((double) pss->count_terminal_segments());
    numtrees++;
  }
  nsd.n = (double) numtrees;
  nsd.mean /= nsd.n; // mean total length of arbor
  termnsd.mean /= termnsd.n; // mean terminal segment length
  if (internsd) internsd->mean /= internsd->n; // mean intermediate segment length
  // [***INCOMPLETE] The std of intermediate segment lengths is not yet computed
  nsd.std = 0.0;
  termnsd.std = 0.0;
  if (numtrees==1) {
    nsd.std = -1.0;
    termnsd.std = -1.0;
  } else {
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1)
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,pss) {
      lendiff = pss->cache - nsd.mean;
      nsd.std += (lendiff*lendiff);
      double termlen;
      PLL_LOOP_FORWARD_NESTED(terminal_segment,pss->TerminalSegments()->head(),1,ts) {
	termlen = ts->Length();
	if (termlen<termnsd.min) termnsd.min = termlen;
	if (termlen>termnsd.max) termnsd.max = termlen;
	lendiff = termlen - termnsd.mean;
	termnsd.std += (lendiff*lendiff);
      }
    }
    // standard SAMPLE deviation sqrt((sumsq(a-mean(a)))/(n-1))
    nsd.std = sqrt(nsd.std/(nsd.n - 1.0));
    termnsd.std = sqrt(termnsd.std/(termnsd.n - 1.0));
  }
  return nsd.mean;
}

double network::mean_number_of_input_terminal_segments(network_statistics_data * nsd) {
  // nt is the mean number of terminal segments per dendritic tree
  // (multipolar neurons have multiple dendritic trees)
  // std can return the standard deviation from that mean
  // arborsamples can return the number of arbors over which the statistic
  // was collected (so that a measure of statistical significance is possible)
  double mean_nt, d_nt, sumsq_nt = 0.0, numtrees = 0.0, nt;
  int sum_nt = 0;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    numtrees += (double) e->InputStructure()->length();
    sum_nt += e->total_input_terminal_segments();
  }
  mean_nt = ((double) sum_nt) / numtrees;
  if (nsd) {
    nsd->n = numtrees; // the number of arbor samples
    nsd->mean = mean_nt;
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1)
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,pss) {
      nt = ((double) (pss->count_terminal_segments()));
      if (nt<nsd->min) nsd->min = nt;
      if (nt>nsd->max) nsd->max = nt;
      d_nt = nt - mean_nt;
      sumsq_nt += (d_nt)*(d_nt);
    }
    if (numtrees==1.0) nsd->std = -1.0; 
    else {
      sumsq_nt /= (numtrees - 1.0);
      nsd->std = sqrt(sumsq_nt); // standard deviation sqrt((sumsq(a-mean(a)))/(n-1))
    }
  }
  return mean_nt;
}

double network::mean_postsynaptic_structure_length(network_statistics_data & termnsd, network_statistics_data & nsd, network_statistics_data * internsd) {
  // This returns the mean length of dendritic arbor in the network and the
  // mean length of dendritic terminal segments in the network.
  double arborlength, totterminalseglength, lendiff;
  int numtrees = 0;
  nsd.mean = 0.0;
  termnsd.n = 0.0;
  termnsd.mean = 0.0;
  if (internsd) {
    internsd->mean = 0.0;
    internsd->n = 0.0;
  }
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1)
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,pss) {
    arborlength = pss->total_arbor_length(&totterminalseglength,internsd);
    if (arborlength<nsd.min) nsd.min = arborlength;
    if (arborlength>nsd.max) nsd.max = arborlength;
    nsd.mean += arborlength;
    termnsd.mean += totterminalseglength;
    termnsd.n += ((double) pss->count_terminal_segments());
    numtrees++;
  }
  nsd.n = (double) numtrees;
  nsd.mean /= nsd.n; // mean total length of arbor
  termnsd.mean /= termnsd.n; // mean terminal segment length
  if (internsd) internsd->mean /= internsd->n; // mean intermediate segment length
  nsd.std = 0.0;
  termnsd.std = 0.0;
  if (numtrees==1) {
    nsd.std = -1.0;
    termnsd.std = -1.0;
  } else {
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1)
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,pss) {
      lendiff = pss->cache - nsd.mean;
      nsd.std += (lendiff*lendiff);
      double termlen;
      PLL_LOOP_FORWARD_NESTED(terminal_segment,pss->TerminalSegments()->head(),1,ts) {
	termlen = ts->Length();
	if (termlen<termnsd.min) termnsd.min = termlen;
	if (termlen>termnsd.max) termnsd.max = termlen;
	lendiff = termlen - termnsd.mean;
	termnsd.std += (lendiff*lendiff);
      }
    }
    // standard SAMPLE deviation sqrt((sumsq(a-mean(a)))/(n-1))
    nsd.std = sqrt(nsd.std/(nsd.n - 1.0));
    termnsd.std = sqrt(termnsd.std/(termnsd.n - 1.0));
  }
  return nsd.mean;
}

void network::move(spatial & pos) {
  // Recenter all neurons on pos.
  // This is useful, for example, as the first step in neurite angle analysis.
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->move(pos);
}

void network::fanin_rot() {
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->fanin_rot();
}

void network::abstract_connections() {
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->abstract_connections();
  progress("Largest abstract connection strength calculated: "+String(max_abstract_strength,"%.3f")+'\n');
}

Fig_Group * network::abstract_connections_Fig() {
  Fig_Group * figgroup = new Fig_Group(netinfo);
  // [***NOTE] If I want to, I can later find an elegant way to delegate
  // much of the following to member functions of the objects involved.
  // I could also move all of the functions that are involved in creating
  // the abstract connections graphs into the connectivity_graph class.
  int j;
  long color, depth, linewidth, linestyle;
  connectivity_graph cg(*this);
  double dashlength = round(cg.graph_R/4000.0), dl;
  if (dashlength<=0.0) dashlength = 1.0;
  for (int i = 0; i<cg.N; i++) {
    color = 0; if (figattr_use_color) color=cg.n[i]->figneuroncolor;
    figgroup->Add_Fig_Object(new Fig_Circle(0,1,color,7,998,-1,0.0,0.0,cg.X[i],cg.Y[i],cg.R[i],cg.X[i],cg.Y[i],cg.X[i]+150,cg.Y[i]+75));
    PLL_LOOP_FORWARD_NESTED(connection,cg.n[i]->OutputConnections()->head(),1,c) if (c->Abstract_Strength()>figattr_connections_threshold) {
      if ((j=cg.find(c->PostSynaptic()))<0) {
	warning("Warning: Postsynaptic neuron of output connection is not in list of neuron in network::abstract_connections_Fig()\n");
      } else {
	if (cg.n[i]->TypeID()==INTERNEURON) color = 4;
	else color = 1;
	double relstrength = c->Abstract_Strength(); // - figattr_connections_threshold;
	depth = i & 511; depth = 997-depth;
	dl = 0.0;
	linestyle = 0;
	linewidth = lround(3.0*fabs(relstrength)/cg.maxstrength);
	if (linewidth<1) {
	  linewidth = 1;
	  linestyle = 1;
	  dl = dashlength;
	}
	figgroup->Add_Fig_Object(new Fig_Line(linestyle,linewidth,color,7,depth,-1,dl,0,0,cg.X[i],cg.Y[i],cg.X[j],cg.Y[j]));
      } 
    }
  }
  return figgroup;
}

Fig_Group * network::Fig_Abstract_Connections(String figname, double width, bool overwrite) {
  // Appends (or overwrites) .fig output to the file figname with a network
  // width specified in micrometers.
  // Note: Regular XFig scale is 1200 dpi. To plot over a width of 6.5
  // inches (1 inch margins on Letter paper) for a network box with specified
  // width in micrometers, that width should equal 6.5*1200 XFig points.
  // E.g. scale 4*292.0 micron = 6.5*1200 => 1 micron = 6.5*1200 / (4*292)
  String netfig;
  if ((overwrite) || (!read_file_into_String(figname,netfig,false))) {
    Fig_Header figheader;
    netfig = figheader.str();
    netfig += colortable->Fig_str();
  }
  if (width>0.0) figscale = (tex_textwidth*1200.0)/width;
  Fig_Group * figgroup = abstract_connections_Fig();
  netfig += figgroup->str();
  write_file_from_String(figname,netfig);
  return figgroup;
}

bool sampling = false; // sampling flag to distinguish callers

Fig_Group * network::net_Fig() {
  Fig_Group * figgroup = new Fig_Group(netinfo);
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) figgroup->Add_Fig_Object(e->net_Fig());
  if ((figattr_show_spatial_segment_subsets) && (sss)) figgroup->Add_Fig_Object(sss->net_Fig());
  switch (figattr_sample_show_progress) {
  case 1: figgroup->Add_Fig_Object(time_Fig(*figgroup)); break;;
  case 2: figgroup->Add_Fig_Object(progress_bar_Fig()); break;;
  }
  if ((figattr_show_scale) && (!sampling)) figgroup->Add_Fig_Object(scale_bar_Fig(*figgroup));
#ifdef VECTOR3D
  if ((figattr_show_axis_arrows) && (!sampling)) figgroup->Add_Fig_Object(axis_arrows_Fig(*figgroup));
#endif
  return figgroup;
}

Fig_Text * network::time_Fig(Fig_Object & rfo) {
  // rfo is a reference object, e.g. a Fig_Group, used to locate the time
  // string Fig_Text
  long days = ((long) t) / 60;
  double seconds = t - ((double) (days*60));
  long minutes = days % 60; days /= 60;
  long hours = days % 24; days /= 24;
  long X = (rfo.TopLeftX()+rfo.BottomRightX())/2, Y=rfo.BottomRightY()+360;
  Fig_Text * timefig = new Fig_Text(1,colortable->colnum(CT_progress_text),1,0,24,0.0,X,Y,String(days)+" days "+String(hours)+" hours "+String(minutes)+" minutes "+String(seconds,"%.3f")+" seconds");
  return timefig;
}

Fig_Group * network::progress_bar_Fig() {
  Fig_Group * proggroup = new Fig_Group();
  double x, y;
#ifdef VECTOR3D
  report("center: "+String(center.X(),"%.3f,")+String(center.Y(),"%.3f,")+String(center.Z(),"%.3f\n"));
#endif
  center.plane_mapped(x,y);
  long unitsy = FIGSCALED(2.7*y);
  long unitsx1 = FIGSCALED(x);
  long unitsx2 = FIGSCALED(2.0*x);
  double days = max_growth_time / 86400.0;
  proggroup->Add_Fig_Object(new Fig_Text(0,colortable->colnum(CT_progress_text),1,0,24,0.0,0,unitsy,"0.0"));
  proggroup->Add_Fig_Object(new Fig_Text(1,colortable->colnum(CT_progress_text),1,0,24,0.0,unitsx1,unitsy,"days"));
  proggroup->Add_Fig_Object(new Fig_Text(2,colortable->colnum(CT_progress_text),1,0,24,0.0,unitsx2,unitsy,String(days,"%.1f")));
  proggroup->Add_Fig_Object(new Fig_Rectangle(0,1,colortable->colnum(CT_progress_text),7,1,-1,0.0,0,0,0,unitsy+300,unitsx2,unitsy+450));
  double progress = t/max_growth_time;
  if (progress>1.0) progress = 1.0;
  long progressx = (long) (progress*((double) unitsx2));
  if ((progressx+30)>=unitsx2) progressx = unitsx2-30;
  proggroup->Add_Fig_Object(new Fig_Rectangle(0,1,colortable->colnum(CT_progress_text),colortable->colnum(CT_progress_text),1,20,0.0,0,0,progressx,unitsy+270,progressx+30,unitsy+480));
  return proggroup;
}

Fig_Group * network::scale_bar_Fig(Fig_Object & rfo) {
  Fig_Group * scalegroup = new Fig_Group();
  long X = rfo.TopLeftX();
  long Y = rfo.BottomRightY()+360;
  long scalebarsize = FIGSCALED(100.0); // a fixed scale bar of 100 micrometers
  scalegroup->Add_Fig_Object(new Fig_Rectangle(0,1,colortable->colnum(CT_progress_text),colortable->colnum(CT_progress_text),1,20,0.0,0,0,X,Y,X+scalebarsize,Y+100));
  scalegroup->Add_Fig_Object(new Fig_Text(0,colortable->colnum(CT_progress_text),1,0,24,0.0,X,Y+400,"100 um"));
  return scalegroup;
}

Fig_Group * network::scale_bar_Fig(long xoffset, long yoffset) {
  Fig_Group * scalegroup = new Fig_Group();
  long X = xoffset;
  long Y = yoffset+360;
  long scalebarsize = FIGSCALED(100.0); // a fixed scale bar of 100 micrometers
  scalegroup->Add_Fig_Object(new Fig_Rectangle(0,1,colortable->colnum(CT_progress_text),colortable->colnum(CT_progress_text),1,20,0.0,0,0,X,Y,X+scalebarsize,Y+100));
  scalegroup->Add_Fig_Object(new Fig_Text(0,colortable->colnum(CT_progress_text),1,0,24,0.0,X,Y+400,"100 um"));
  return scalegroup;
}

Fig_Group * network::axis_arrows_Fig(Fig_Object & rfo) {
  Fig_Group * axisarrowsgroup = new Fig_Group();
  long X1, Y1, X2, Y2;
  long Xoffset = rfo.TopLeftX();
  long Yoffset = rfo.TopLeftY();
  //spatial minextent, maxextent;
  //spatial_extents(minextent,maxextent);
  double axisscale = 0.1 * INV_FIGSCALED((rfo.BottomRightX()-rfo.TopLeftX()));
  spatial origin(0.0,0.0,0.0);
  spatial xaxis(axisscale,0.0,0.0);
  spatial yaxis(0.0,axisscale,0.0);
  spatial zaxis(0.0,0.0,axisscale);
  double x,y;
  origin.plane_mapped(x,y);
  X1 = FIGSCALED(x);
  Y1 = FIGSCALED(y);
  xaxis.plane_mapped(x,y);
  X2 = FIGSCALED(x);
  Y2 = FIGSCALED(y);
  axisarrowsgroup->Add_Fig_Object(new Fig_Line(0,2,colortable->colnum(CT_progress_text),colortable->colnum(CT_progress_text),1,-1,0.0,0,0,Xoffset+X1,Yoffset+Y1,Xoffset+X2,Yoffset+Y2));
  axisarrowsgroup->Add_Fig_Object(new Fig_Text(0,colortable->colnum(CT_progress_text),1,0,24,0.0,Xoffset+X2+15,Yoffset+Y2,"x"));
  yaxis.plane_mapped(x,y);
  X2 = FIGSCALED(x);
  Y2 = FIGSCALED(y);
  axisarrowsgroup->Add_Fig_Object(new Fig_Line(0,2,colortable->colnum(CT_progress_text),colortable->colnum(CT_progress_text),1,-1,0.0,0,0,Xoffset+X1,Yoffset+Y1,Xoffset+X2,Yoffset+Y2));
  axisarrowsgroup->Add_Fig_Object(new Fig_Text(0,colortable->colnum(CT_progress_text),1,0,24,0.0,Xoffset+X2+15,Yoffset+Y2,"y"));
  zaxis.plane_mapped(x,y);
  X2 = FIGSCALED(x);
  Y2 = FIGSCALED(y);
  axisarrowsgroup->Add_Fig_Object(new Fig_Line(0,2,colortable->colnum(CT_progress_text),colortable->colnum(CT_progress_text),1,-1,0.0,0,0,Xoffset+X1,Yoffset+Y1,Xoffset+X2,Yoffset+Y2));
  axisarrowsgroup->Add_Fig_Object(new Fig_Text(0,colortable->colnum(CT_progress_text),1,0,24,0.0,Xoffset+X2+15,Yoffset+Y2,"z"));
  return axisarrowsgroup;
}

Fig_Group * network::Fig_Output(String figname, double width, bool overwrite) {
  // Appends (or overwrites) .fig output to the file figname with a network
  // width specified in micrometers.
  // Note: Regular XFig scale is 1200 dpi. To plot over a width of 6.5
  // inches (1 inch margins on Letter paper) for a network box with specified
  // width in micrometers, that width should equal 6.5*1200 XFig points.
  // E.g. scale 4*292.0 micron = 6.5*1200 => 1 micron = 6.5*1200 / (4*292)
  String netfig;
  XFIG_RANGE_CHECK_INIT();
  if ((overwrite) || (!read_file_into_String(figname,netfig,false))) {
    Fig_Header figheader;
    netfig = figheader.str();
    netfig += colortable->Fig_str();
  }
  if (width>0.0) figscale = (tex_textwidth*1200.0)/width;
  if (!netinfo.contains("figscale")) {
    netinfo += "\nfigscale=" + String(figscale,"%.3f\n");
    netinfo += "Normal XFig scale is 1200 dpi\n";
  }
  Fig_Group * figgroup = net_Fig();
  netfig += figgroup->str();
  write_file_from_String(figname,netfig);
  if (Fig_output_modified>0) warning(XFIG_RANGE_CHECK_WARNING(figname,tex_textwidth));
  return figgroup;
}

Fig_Group * network::Fig_Neurons(String figname, double width) {
  // Overwrites .fig output to a set of files with the base name
  // figname, in which XXXXX is replaced with neuron indices, using
  // a network width specified in micrometers (see Fig_Output()).
  Fig_Header figheader;
  XFIG_RANGE_CHECK_INIT();
  if (width>0.0) figscale = (tex_textwidth*1200.0)/width;
  if (!netinfo.contains("figscale")) netinfo += "figscale=" + String(figscale,"%.3f\n");
  long cnt = 0;
  // [***NOTE] Efficiency gains can be made if this function is called before Fig_Output(),
  // and if this function also places the fig_objects created with e->net_Fig() calls into
  // a combined Fig_Group. Then Fig_Output() can test if that Fig_Group is already provided,
  // instead of having to create all the objects again.
  Fig_Group * figgroup = NULL;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    String netfig(figheader.str());
    netfig += colortable->Fig_str();
    if (figgroup) delete figgroup;
    figgroup = new Fig_Group(netinfo);
    figgroup->Add_Fig_Object(e->net_Fig());
    netfig += figgroup->str();
    String indexedfigname(figname); String cntstr(cnt); while (cntstr.length()<5) cntstr.prepend('0');
    indexedfigname.gsub("XXXXX",cntstr);
    write_file_from_String(indexedfigname,netfig);
    cnt++;
  }
  if (Fig_output_modified>0) warning(XFIG_RANGE_CHECK_WARNING(figname,tex_textwidth));
  return figgroup;
}

Txt_Group * network::net_Txt() {
  // Initialize strings and pre-allocate space
  if (outattr_Txt_separate_files) {
    if (NES_output) {
      Txt_neuronlist.reset(new String(COLUMN_LABELS_NEURONS_WITH_RADIUS));
    } else {
      Txt_neuronlist.reset(new String(COLUMN_LABELS_NEURONS));
    }
    Txt_synapselist.reset(new String(COLUMN_LABELS_SYNAPSES));
  } else {
    Txt_neuronlist.reset(new String());
    Txt_synapselist.reset(new String());
  }
  Txt_neuronindex = 0;
  Txt_synapseindex = 0;
  int nlsize = (PLLRoot<neuron>::length()+2)*64; // estimated from measurements in output
  Txt_neuronlist->alloc(nlsize);
  long unsigned int totnumsynapses = 0;
  for (synapse_type i = synapse_type(0); i<syntype_IDs; i=synapse_type(i+1)) totnumsynapses += (synapse_inventory[i][SYNGENESIS] - synapse_inventory[i][SYNLOSS]);
  Txt_synapselist->alloc((totnumsynapses+3)*128); // estimated from measurements in output
  if (outattr_track_nodegenesis) {
    Txt_nodeindex = -1;
    if (outattr_Txt_separate_files) {
      Txt_fiberrootlist.reset(new String(COLUMN_LABELS_ROOT_NODES));
      Txt_continuationnodelist.reset(new String(COLUMN_LABELS_FIBER_NODES));
      Txt_bifurcationnodelist.reset(new String(COLUMN_LABELS_FIBER_NODES));
      Txt_terminalgrowthconelist.reset(new String(COLUMN_LABELS_FIBER_NODES));
    } else {
      Txt_fiberrootlist.reset(new String());
      Txt_continuationnodelist.reset(new String());
      Txt_bifurcationnodelist.reset(new String());
      Txt_terminalgrowthconelist.reset(new String());
    }
    long unsigned int totalnumroots = 0, continuation = 0, bifurcation = 0, terminal = 0;
    PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
      totalnumroots += e->InputStructure()->length() + e->OutputStructure()->length();
      parsing_fs_type = dendrite_fs;
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,d) d->count_segment_types(continuation,bifurcation,terminal);
      parsing_fs_type = axon_fs;
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,a) a->count_segment_types(continuation,bifurcation,terminal);
    }
    Txt_fiberrootlist->alloc((totalnumroots+3)*80);
    Txt_continuationnodelist->alloc((continuation+3)*128);
    Txt_bifurcationnodelist->alloc((bifurcation+3)*128);
    Txt_terminalgrowthconelist->alloc((terminal+3)*128);
  }

  // Initialize a text output group
  Txt_Group * txtgroup = new Txt_Group(netinfo);

  // Collect information from all neurons in the network
  //***  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->net_Txt();
  progress("Creating structure output files in text format:\n");
  int nnum = 0;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    progress("neuron "+String((long) nnum));
    e->net_Txt();
    nnum++;
  }
  //PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) txtgroup->Add_Txt_Object(e->net_Txt());
  return txtgroup;
}

bool network::Txt_Output(String txtname) {
  // Overwrites .txt output to the file figname.
  // This can be used as a simple function to obtain the structure data of the network, not its state at some time during a simulation.
  if (outattr_Txt_separate_files) {
    Txt_Header txtheader;
    String fullheader(txtheader.str()+net_Txt()->str());
    write_file_from_String(txtname+".header",fullheader);
    write_file_from_String(txtname+".neurons",*(Txt_neuronlist));
    //Txt_neuronlist.reset();
    write_file_from_String(txtname+".synapses",*(Txt_synapselist));
    //Txt_synapselist.reset();
    if (outattr_track_nodegenesis) {
      write_file_from_String(txtname+".rootnodes",*(Txt_fiberrootlist));
      //Txt_fiberrootlist.reset();
      write_file_from_String(txtname+".continuationnodes",*(Txt_continuationnodelist));
      //Txt_continuationnodelist.reset();
      write_file_from_String(txtname+".bifurcationnodes",*(Txt_bifurcationnodelist));
      //Txt_bifurcationnodelist.reset();
      write_file_from_String(txtname+".growthcones",*(Txt_terminalgrowthconelist));
      //Txt_terminalgrowthconelist.reset();
      write_file_from_String(txtname+".tuftnodes",*(Txt_tuftrootbranchnodelist));
      write_file_from_String(txtname+".obliquenodes",*(Txt_obliquerootbranchnodelist));
    }
  } else {
    String nettxt;
    Txt_Header txtheader;
    //  if (!read_file_into_String(txtname,nettxt,false)) {
    nettxt = txtheader.str();
    //   }
    nettxt += net_Txt()->str(); // processed 'netinfo'
    nettxt += "neurons:\n";
    nettxt += (*Txt_neuronlist);
    //Txt_neuronlist.reset();
    nettxt += "synapses:\n";
    nettxt += (*Txt_synapselist);
    //Txt_synapselist.reset();
    if (outattr_track_nodegenesis) {
      nettxt += "fiber structure root nodes:\n";
      nettxt += (*Txt_fiberrootlist);
      //Txt_fiberrootlist.reset();
      nettxt += "fiber continuation nodes:\n";
      nettxt += (*Txt_continuationnodelist);
      //Txt_continuationnodelist.reset();
      nettxt += "fiber bifurcation nodes:\n";
      nettxt += (*Txt_bifurcationnodelist);
      //Txt_bifurcationnodelist.reset();
      nettxt += "terminal fiber growth cones:\n";
      nettxt += (*Txt_terminalgrowthconelist);
      //Txt_terminalgrowthconelist.reset();
      nettxt += "apical dendrite tuft root nodes:\n";
      nettxt += (*Txt_tuftrootbranchnodelist);
      nettxt += "apical dendrite oblique root nodes:\n";
      nettxt += (*Txt_obliquerootbranchnodelist);
    }
    write_file_from_String(txtname,nettxt);
  }
  return true;
}

/**
 * Generate output data structures used to deliver the grown model
 * to NES (or another platform) within which Netmorph is embedded.
 */
bool network::net_NES() {
  // *** NOTE: At present, including only those parts used in NES model building
  //     as drafted in the Python frontend scripts netmorph2nes.py and netmorph2obj.py.
  //     These use collected data about:
  //     - neurons
  //     - synapses
  //     - bifurcationnodes
  //     - continuationnodes
  //     - growthcones
  //     - rootnodes
  std::unique_ptr<NetData> netdata = std::make_unique<NetData>();
  progress("Collecting structure output data in NES format:\n");
  int nnum = 0;
  PLL_LOOP_FORWARD(neuron, PLLRoot<neuron>::head(), 1) {
    progress("neuron "+String((long) nnum));
    if (!e->net_NES(*netdata)) {
      return false;
    }
    nnum++;
  }
  return true;
}

bool network::Data_Output_Synapse_Distance(String dataname) {
  int numbins = (int) (10000.0 / outattr_distance_frequency_distbinsize);
  unsigned long * axonbins = new unsigned long[numbins];
  unsigned long * dendritebins =  new unsigned long[numbins];
  for (int i = 0; i<numbins; i++) axonbins[i] = dendritebins[i] = 0;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    PLL_LOOP_FORWARD_NESTED(connection,e->OutputConnections()->head(),1,c) {
      PLL_LOOP_FORWARD_NESTED(synapse,c->Synapses()->head(),1,s) {
        int axonbinnum = (int) (c->PreSynaptic()->Pos().distance(s->Structure()->P0) / outattr_distance_frequency_distbinsize);
        int dendritebinnum = (int) (c->PostSynaptic()->Pos().distance(s->Structure()->P0) / outattr_distance_frequency_distbinsize);
	if ((axonbinnum>=numbins) || (dendritebinnum>=numbins)) {
	  int maxbinnum = axonbinnum;
	  if (dendritebinnum>maxbinnum) maxbinnum = dendritebinnum;
	  int newnumbins = numbins + 2*((maxbinnum-numbins)+1);
	  unsigned long * oldaxonbins = axonbins;
	  unsigned long * olddendritebins = dendritebins;
	  axonbins = new unsigned long[newnumbins];
	  dendritebins = new unsigned long[newnumbins];
	  for (int i = 0; i<numbins; i++) {
	    axonbins[i] = oldaxonbins[i];
	    dendritebins[i] = olddendritebins[i];
	  }
	  for (int i = numbins; i<newnumbins; i++) axonbins[i] = dendritebins[i] = 0;
	  numbins = newnumbins;
	  delete[] oldaxonbins;
	  delete[] olddendritebins;
	}
	axonbins[axonbinnum]++;
	dendritebins[dendritebinnum]++;
      }
    }
  }
  String datastr("% "+dataname+"\n% NETMORPH output: Synapses by distance (in micrometers)\n\n");
  datastr += "distbinsize = " + String(outattr_distance_frequency_distbinsize,"%.3f;\n");
  datastr += "numbins = " + String((long) numbins) + ";\n";
  datastr += "% Radial distance from presynaptic neuron to synapse location\naxonbins = [\n";
  for (int i = 0; i<numbins; i++) datastr += String((long) axonbins[i]) + '\n';
  datastr += "];\n";
  datastr += "% Radial distance from postsynaptic neuron to synapse location\ndendritebins = [\n";
  for (int i = 0; i<numbins; i++) datastr += String((long) dendritebins[i]) + '\n';
  datastr += "];\n";
  write_file_from_String(dataname,datastr);
  delete[] axonbins;
  delete[] dendritebins;
  return true;
}

bool network::Data_Output_Connection_Distance(String dataname) {
  int numbins = (int) (10000.0 / outattr_distance_frequency_distbinsize);
  unsigned long * bins = new unsigned long[numbins];
  for (int i = 0; i<numbins; i++) bins[i] = 0;
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) {
    PLL_LOOP_FORWARD_NESTED(connection,e->OutputConnections()->head(),1,c) {
      int binnum = (int) (c->PreSynaptic()->Pos().distance(c->PostSynaptic()->Pos()) / outattr_distance_frequency_distbinsize);
      if (binnum>=numbins) {
	int newnumbins = numbins + 2*((binnum-numbins)+1);
	unsigned long * oldbins = bins;
	bins = new unsigned long[newnumbins];
	for (int i = 0; i<numbins; i++) bins[i] = oldbins[i];
	for (int i = numbins; i<newnumbins; i++) bins[i] = 0;
	numbins = newnumbins;
	delete[] oldbins;
      }
      bins[binnum]++;
    }
  }
  String datastr("% "+dataname+"\n% NETMORPH output: Connections by distance (in micrometers)\n\n");
  datastr += "distbinsize = " + String(outattr_distance_frequency_distbinsize,"%.3f;\n");
  datastr += "numbins = " + String((long) numbins) + ";\n";
  datastr += "bins = [\n";
  for (int i = 0; i<numbins; i++) datastr += String((long) bins[i]) + '\n';
  datastr += "];\n";
  write_file_from_String(dataname,datastr);
  delete[] bins;
  return true;
}

VRML_Group * network::net_VRML() {
  if (VRML_neuronlist) delete VRML_neuronlist;
  if (VRML_synapselist) delete VRML_synapselist;
  VRML_neuronlist = new String();
  VRML_synapselist = new String();
  VRML_neuronindex = 0;
  VRML_synapseindex = 0;
  int nlsize = PLLRoot<neuron>::length()*64;
  VRML_neuronlist->alloc(nlsize);
  VRML_synapselist->alloc(nlsize*64); // pre-allocate some space to save time
  VRML_nodeindex = -1;
  if (VRML_fiberlist) delete VRML_fiberlist;
  VRML_fiberlist = new String();
  VRML_fiberlist->alloc(nlsize*8*64);
  if (VRML_continuationnodelist) delete VRML_continuationnodelist;
  VRML_continuationnodelist = new String();
  VRML_continuationnodelist->alloc(nlsize*64*64);
  if (VRML_bifurcationnodelist) delete VRML_bifurcationnodelist;
  VRML_bifurcationnodelist = new String();
  VRML_bifurcationnodelist->alloc(nlsize*32*64);
  if (VRML_terminalgrowthconelist) delete VRML_terminalgrowthconelist;
  VRML_terminalgrowthconelist = new String();
  VRML_terminalgrowthconelist->alloc(nlsize*32*64);
  VRML_Group * vrmlgroup = new VRML_Group(netinfo);
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->net_VRML();
  //PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) txtgroup->Add_VRML_Object(e->net_VRML());
  return vrmlgroup;
}

bool network::VRML_Output(String vrmlname) {
  // Overwrites .x3d output to the file vrmlname.
  // This is intended to produce VRML compatible 3D output. From there,
  // conversion to other file formats can take place.
  String netvrml;
  VRML_Header vrmlheader;
  netvrml = vrmlheader.str();
  netvrml += net_VRML()->str(); // processed 'netinfo'
  netvrml += (*VRML_neuronlist);
  delete VRML_neuronlist; VRML_neuronlist = NULL;
  netvrml += (*VRML_synapselist);
  delete VRML_synapselist; VRML_synapselist = NULL;
  netvrml += (*VRML_fiberlist);
  delete VRML_fiberlist; VRML_fiberlist = NULL;
  netvrml += (*VRML_continuationnodelist);
  delete VRML_continuationnodelist; VRML_continuationnodelist = NULL;
  netvrml += (*VRML_bifurcationnodelist);
  delete VRML_bifurcationnodelist; VRML_bifurcationnodelist = NULL;
  netvrml += (*VRML_terminalgrowthconelist);
  delete VRML_terminalgrowthconelist; VRML_terminalgrowthconelist = NULL;
  //netvrml += "<NavigationInfo type='\"EXAMINE\"'/>\n</Scene>\n</X3D>\n";
  netvrml += "</Scene>\n</X3D>\n";
  write_file_from_String(vrmlname,netvrml);
  return true;
}

Catacomb_projectionsmatrix * cpm = NULL;

Catacomb_Group * network::net_Catacomb() {
  // [***NOTE] This function may allocate a Catacomb_projectionsmatrix. That
  // can be unallocated in the calling function.
  Catacomb_neuronindex = 0;
  Catacomb_connectionindex = 0;
  Catacomb_Group * ccmgroup = new Catacomb_Group(netinfo);
  // 1. Collect information about regions (one or more neurons)
  if (ccmregionsgroup) delete ccmregionsgroup;
  ccmregionsgroup = new Catacomb_Group();
  PLL_LOOP_FORWARD(region,regions.head(),1) e->net_Catacomb();
  // 2. Collect information about projections and connection tables within and between regions
  if (ccmprojectionsgroup) delete ccmprojectionsgroup;
  ccmprojectionsgroup = new Catacomb_Group();
  if (ccmconnectionsgroup) delete ccmconnectionsgroup;
  ccmconnectionsgroup = new Catacomb_Group();
  // create a matrix of regions to regions
  if (regions.length()>0) {
    cpm = new Catacomb_projectionsmatrix(regions.length());
    int i=0;
    PLL_LOOP_FORWARD(region,regions.head(),1) { // in each region
      PLL_LOOP_FORWARD_NESTED(neuronptrlist,e->Nlist().head(),1,nptrlist) { // for each neuron
	PLL_LOOP_FORWARD_NESTED(connection,nptrlist->N()->OutputConnections()->head(),1,c) { // check all connections
	  neuron * n_post = c->PostSynaptic();
	  int j=regions.find(*n_post);
	  if (j>=0) {
	    if (cpm->el(i,j)==NULL) cpm->el(i,j) = new Catacomb_connectiontable(e->Nlist().length(),regions[j]->Nlist().length());
	    cpm->el(i,j)->el(e->find(*nptrlist->N()),regions[j]->find(*n_post)) = 1;
	  }
	}
      }
      i++;
    }
    ccmprojectionsgroup->Add_Catacomb_Object(new Catacomb_Projection(*cpm,regions)); // [***NOTE] Alternatively, one object can be created for each projection or connection table.
    ccmconnectionsgroup->Add_Catacomb_Object(new Catacomb_ConnectionTable(*cpm,regions));
  }
  // link all the data structures in the desired order
  ccmgroup->Add_Catacomb_Object(ccmregionsgroup);
  ccmgroup->Add_Catacomb_Object(ccmprojectionsgroup);
  Catacomb_Neuron * ccmneurondef = new Catacomb_Neuron();
  ccmgroup->Add_Catacomb_Object(ccmneurondef);
  ccmgroup->Add_Catacomb_Object(ccmconnectionsgroup);
  return ccmgroup;
}

bool network::Catacomb_Output(String ccmname) {
  // Overwrites .ccm output to the file ccmname.
  // This is intended to produce Catacomb compatible output.
  // For the .ccm format, produce output as follows:
  //   header (includes fixedsteprunner)
  //     list of network + vectorrecorder + joinclonerelay + spike cables + vector cables
  //     list of spike projections
  //     close assembly items and workbench and list
  //     neuron definition
  //     define connection tables
  //   close list and file
  String netccm;
  Catacomb_Header ccmheader;
  netccm = ccmheader.str();
  Catacomb_Group * ccmnet = net_Catacomb();
  //cout << "SEGFAULT SEARCH: @5\n"; cout.flush();
  netccm += ccmnet->str(); // processed 'netinfo'
  //cout << "SEGFAULT SEARCH: @6\n"; cout.flush();
  if (cpm) delete cpm;
  delete ccmnet;
  netccm += "</lists>\n<name>\"box\"</name>\n</org.enorg.catacomb.core.CcmbBox>\n";
  write_file_from_String(ccmname,netccm);
  return true;
}

#ifdef VECTOR3D
void network::net_Slice(Slice & slice) { // [***INCOMPLETE] Change this to a different return type.
  // The iterator functions of the Slice class are used to walk through all defined slices and return those that are single slices with vertices, which are then given to the neuron::net_Slice() function. Iteration is done *within* the slice objects.
  for (Slice * s = slice.iterator_first(); (s); s = s->iterator_next()) PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->net_Slice(s);
}

bool network::Slice_Output(String slicename, Slice & slice) {
  // Generates output files that simulate the appearance of histological slices.
  // [***INCOMPLETE] I need to decide which file format to use, then store the output of net_Slice().
  net_Slice(slice);
  return true;
}

Fig_Group * network::Slice_Outlines_Output(String figname, Slice & slice, double width, bool overwrite) {
  // Appends (or overwrites) .fig output to the file figname with a network
  // width specified in micrometers.
  // Note: Regular XFig scale is 1200 dpi. To plot over a width of 6.5
  // inches (1 inch margins on Letter paper) for a network box with specified
  // width in micrometers, that width should equal 6.5*1200 XFig points.
  // E.g. scale 4*292.0 micron = 6.5*1200 => 1 micron = 6.5*1200 / (4*292)
  String netfig;
  XFIG_RANGE_CHECK_INIT();
  if ((overwrite) || (!read_file_into_String(figname,netfig,false))) {
    Fig_Header figheader;
    netfig = figheader.str();
    netfig += colortable->Fig_str();
  }
  if (width>0.0) figscale = (tex_textwidth*1200.0)/width;
  if (!netinfo.contains("figscale")) {
    netinfo += "\nfigscale=" + String(figscale,"%.3f\n");
    netinfo += "Normal XFig scale is 1200 dpi\n";
  }
  Fig_Group * figgroup = slice.Outlines_Fig();
  netfig += figgroup->str();
  write_file_from_String(figname,netfig);
  if (Fig_output_modified>0) warning(XFIG_RANGE_CHECK_WARNING(figname,tex_textwidth));
  return figgroup;
}
// VECTOR3D
#endif

void network::set_figattr(int fattr) {
  // set the figattr parameter for all neurons in the network
  PLL_LOOP_FORWARD(neuron,PLLRoot<neuron>::head(),1) e->figattr = fattr;
}

void network::visible_pre_and_target_post(neuron & n) {
  // make the presynaptic connections of n visible and the postsynaptic
  // connections of all targets of n, as well as the cell bodies of n and
  // its output targets
  n.figattr = (n.figattr & FIGATTR_SHOWPRESYN & FIGATTR_SHOWCELL) | FIGATTR_REPORTDEPTH;
  PLL_LOOP_FORWARD(connection,n.OutputConnections()->head(),1) e->PostSynaptic()->figattr &= (FIGATTR_SHOWPOSTSYN & FIGATTR_SHOWCELL);
}

/**
 * For all neurons in the network, traverse the tree of fibers
 * down to the terminal segments and at each fiber segment apply
 * the operation specified by the fibre_tree_op object.
 */
void network::tree_op(fibre_tree_op& op) {
  PLL_LOOP_FORWARD(neuron, PLLRoot<neuron>::head(), 1) e->tree_op(op);
}

/**
 * For all neurons in the network, apply the operation specified
 * by the neuron_list_op object.
 */
void network::neuron_op(neuron_list_op& op) {
  PLL_LOOP_FORWARD(neuron, PLLRoot<neuron>::head(), 1) e->neuron_op(op);
}

/**
 * For all neurons in the network, traverse the set of connections
 * and subset of synapses of each connection, and apply the
 * operation specified by the synapse_tree_op object.
 */
void network::synapse_op(synapse_tree_op& op) {
  PLL_LOOP_FORWARD(neuron, PLLRoot<neuron>::head(), 1) e->synapse_op(op);
}

const String cube_header = R"HEAD(# netmorph_output.obj
#

)HEAD";

const String cube_top_begin = R"TOP(o soma)TOP";
const String cube_top_end = R"TOP(
mtllib cube.mtl
)TOP";


const String dendrites_top_begin = R"TOP(o dendrites)TOP";
const String dendrites_top_end = R"TOP(
)TOP";

const String axon_top_begin = R"TOP(o axon)TOP";
const String axon_top_end = R"TOP(
)TOP";

const String faces = R"FACES(
vn 0.000000 0.000000 1.000000
vn 0.000000 1.000000 0.000000
vn 0.000000 0.000000 -1.000000
vn 0.000000 -1.000000 0.000000
vn 1.000000 0.000000 0.000000
vn -1.000000 0.000000 0.000000
)FACES";

const String simple_cube_faces = R"FACES(
vn 0.000000 -1.000000 0.000000
vn 0.000000 1.000000 0.000000
vn 1.000000 0.000000 0.000000
vn -0.000000 0.000000 1.000000
vn -1.000000 -0.000000 -0.000000
vn 0.000000 0.000000 -1.000000
)FACES";

const String cube_extras = R"EXTRAS(
usemtl cube
s off
)EXTRAS";

const int face_lines[6][4][2] = {
  {{0,0}, {1,0}, {2,0}, {3,0}},
  {{4,1}, {7,1}, {6,1}, {5,1}},
  {{0,2}, {4,2}, {5,2}, {1,2}},
  {{1,3}, {5,3}, {6,3}, {2,3}},
  {{2,4}, {6,4}, {7,4}, {3,4}},
  {{4,5}, {0,5}, {3,5}, {7,5}}
};

#define PYRAMID_FACES 6
const int pyramid_face_lines[PYRAMID_FACES][3] = {
  {3, 0, 1},
  {2, 3, 1},
  {4, 1, 0},
  {3, 4, 0},
  {2, 4, 3},
  {4, 2, 1}
};

#define TWO_SIDED_CUBE_FACES 12
const int two_sided_cube_faces[TWO_SIDED_CUBE_FACES][3] = {
  {1, 2, 3},
  {7, 6, 5},
  {4, 5, 1},
  {5, 6, 2},
  {2, 6, 7},
  {0, 3, 7},
  {0, 1, 3},
  {4, 7, 5},
  {0, 4, 1},
  {1, 5, 2},
  {3, 2, 7},
  {4, 0, 7}
};

const std::vector<spatial> simple_cube_lines = {
  { 1.0, -1.0, -1.0 },
  { 1.0, -1.0, 1.0 },
  { -1.0, -1.0, 1.0 },
  { -1.0, -1.0, -1.0 },
  { 1.0, 1.0, -1.0 },
  { 1.0, 1.0, 1.0 },
  { -1.0, 1.0, 1.0 },
  { -1.0, 1.0, -1.0 }
};

const std::vector<spatial> pyramid = {
  {0.0, 0.0, 0.0},
  {1.0, 0.0, 0.0},
  {1.0, 1.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.5, 0.5, 1.6}
};

String obj_vertex(const spatial& v, const String& extra = "") {
  String objstr("v");
  objstr += String(v.X(), " %.3f");
  objstr += String(v.Y(), " %.3f");
  objstr += String(v.Z(), " %.3f");
  objstr += extra;
  objstr += '\n';
  return objstr;
}

String obj_faces(int start, const int* face_id) {
  String objstr("f");
  for (int j = 0; j < 3; j++) {
    objstr += ' '+String(long(start+face_id[j]));
  }
  objstr += '\n';
  return objstr;
}

String obj_line(long v1_idx, long v2_idx) {
  String objstr("l");
  objstr += ' '+String(v1_idx);
  objstr += ' '+String(v2_idx);
  objstr += '\n';
  return objstr;
}

struct ObjData {
  network& net;
  spatial center;
  std::vector<String> soma_obj_names;
  std::vector<String> axon_obj_names;
  std::vector<String> dendrite_obj_names;
  std::vector<String> synapse_obj_names;
  int cube_num = 0;
  int synapse_idx = 0;
  int face_start = 1;
  int vertex_start = 1;
  ObjData(network* _net): net(*_net) {}
};

String add_soma_pyramid(ObjData& data, neuron* n) {
  String pyramid_data(cube_top_begin + String(long(data.cube_num)) + cube_top_end);

  double r = n->Radius();
  double scale = 2.0*r;
  spatial offset = n->Pos() - spatial(r, r, r) - data.center;
  for (auto& line : pyramid) {
    spatial vertex = scale*line + offset;
    pyramid_data += obj_vertex(vertex);
  }
  pyramid_data += faces;
  pyramid_data += cube_extras;

  for (size_t i = 0; i < PYRAMID_FACES; i++) {
    pyramid_data += obj_faces(data.vertex_start, pyramid_face_lines[i]); // face_start or vertex_start?
  }

  data.vertex_start += pyramid.size();
  data.face_start += PYRAMID_FACES;

  return pyramid_data;
}

String add_neuron_neurites(ObjData& data, const String& somalabel, bool doaxon, neuron* n) {
  String neurite_segments;

  class Neurites2Obj: public fibre_tree_op {
  public:
    ObjData& data;
    String& neurite_segments;
    Neurites2Obj(ObjData& _data, String& _neuritesegments): data(_data), neurite_segments(_neuritesegments) {}
    // Selection of neuron an pre- or postsynaptic fiber structure happens
    // outside of this class. See below.
    virtual void op(fibre_segment* fs) {
      spatial v1 = fs->P0 - data.center;
      spatial v2 = fs->P1 - data.center;
      neurite_segments += obj_vertex(v1, " 1.0 0.0 0.0");
      neurite_segments += obj_vertex(v2, " 1.0 0.0 0.0");
      neurite_segments += obj_line(data.vertex_start, data.vertex_start+1);
      data.vertex_start += 2;
    }
  };

  Neurites2Obj neurites2obj(data, neurite_segments);
  neurites2obj.do_dendrites = !doaxon;
  neurites2obj.do_axons = doaxon;
  n->tree_op(neurites2obj);

  return neurite_segments;
}

String add_simple_cube(ObjData& data, const String& synapsestr, double radius, spatial& receptor_loc) {
  String cube_data("o "+synapsestr+'\n');

  spatial offset = receptor_loc - data.center;
  for (auto& line : simple_cube_lines) {
    spatial vertex = radius*line + offset;
    cube_data += obj_vertex(vertex);
  }

  for (size_t i = 0; i < TWO_SIDED_CUBE_FACES; i++) {
    cube_data += obj_faces(data.vertex_start, two_sided_cube_faces[i]); // face_start or vertex_start?
  }

  data.vertex_start += simple_cube_lines.size();
  data.face_start += TWO_SIDED_CUBE_FACES;

  return cube_data;
}

String make_Wavefront_OBJ(ObjData& data) {
  String output(cube_header);

  class Somas2Obj: public neuron_list_op {
  public:
    ObjData& data;
    String& output;
    Somas2Obj(ObjData& _data, String& _output): data(_data), output(_output) {}
    virtual void op(neuron* n) {
      String cubenumstr(std::to_string(data.cube_num).c_str());

      data.soma_obj_names.emplace_back("soma"+cubenumstr);
      output += add_soma_pyramid(data, n);

      data.axon_obj_names.emplace_back("axon"+cubenumstr);
      output += axon_top_begin + cubenumstr + axon_top_end;
      output += '\n' + add_neuron_neurites(data, n->label(), true, n);

      data.dendrite_obj_names.emplace_back("dendrites"+cubenumstr);
      output += dendrites_top_begin + cubenumstr + dendrites_top_end;
      output += '\n' + add_neuron_neurites(data, n->label(), false, n);

      data.cube_num++;
    }
  };

  Somas2Obj somas2obj(data, output);
  data.net.neuron_op(somas2obj);

  class Synapse2Obj: public synapse_tree_op {
  public:
    ObjData& data;
    String& output;
    Synapse2Obj(ObjData& _data, String& _output): data(_data), output(_output) {}
    virtual void op(synapse* s) {
        String synapsestr(("synapse"+std::to_string(data.synapse_idx)).c_str());

        spatial receptor_loc = s->Structure()->P1;
        spatial spine_loc = s->Structure()->P0;

        double synapse_size = norm(receptor_loc - spine_loc);

        data.synapse_obj_names.emplace_back(synapsestr);

        output += '\n' + add_simple_cube(data, synapsestr, synapse_size/2.0, receptor_loc);

        data.synapse_idx++;
    }
  };

  Synapse2Obj synapse2obj(data, output);
  data.net.synapse_op(synapse2obj);

  return output;
}

class FindCenter: public neuron_list_op {
public:
  spatial center;
  size_t num = 0;
  virtual void op(neuron* n) {
    center += n->Pos();
    num++;
  }
  spatial get_center() {
    center /= double(num);
    return center;
  }
};

String VecOfStrings2JSONList(const std::vector<String>& vec) {
  String jsonlist("[");
  for (const auto& s : vec) {
    jsonlist += '"'+s+"\",";
  }
  if (vec.empty()) jsonlist += ']';
  else jsonlist[jsonlist.length()-1] = ']';
  return jsonlist;
}

/**
 * Generate Wavefront OBJ output file containing the network
 * morphology.
 */
bool network::Obj_Output(String objpath, double axonbevdepth, double dendritebevdepth) {

  ObjData data(this);

  FindCenter find_center;
  neuron_op(find_center);
  data.center = find_center.get_center();

  String obj_data = make_Wavefront_OBJ(data);

  write_file_from_String(objpath, obj_data);

  String blender_obj_data = "{\n"
    "\t\"obj_path\": \""+objpath+"\",\n"
    "\t\"axons\": "+VecOfStrings2JSONList(data.axon_obj_names)+",\n"
    "\t\"dendrites\": "+VecOfStrings2JSONList(data.dendrite_obj_names)+",\n"
    "\t\"somas\": "+VecOfStrings2JSONList(data.soma_obj_names)+",\n"
    "\t\"synapses\": "+VecOfStrings2JSONList(data.synapse_obj_names)+" ,\n"
    "\t\"blend_path\": \""+objpath+".blend\",\n"
    "\t\"axonbevdepth\": "+String(axonbevdepth, "%.3f")+",\n"
    "\t\"dendritebevdepth\": "+String(dendritebevdepth, "%.3f")+"\n"
    "}\n";

  write_file_from_String(objpath+".json", blender_obj_data);
  return true;
}

bool network::Blend_Output(String objpath) {
  // 1. Load template Blender Python file and create specific temporary one.
  String objcurves_content(objcurves_py);
  // if (!read_file_into_String("../Source/objcurves.py", objcurves_content)) {
  //   warning("Unable to read Blender Python template file objcurves.py.");
  //   return false;
  // }
  int objdata_pos = objcurves_content.index("OBJDATA_JSON");
  if (objdata_pos<0) {
    warning("OBJDATA_JSON not found in objcurves.py.");
    return false;
  }
  String specific_objcurves(objcurves_content.before(objdata_pos));
  specific_objcurves += objpath+".json";
  specific_objcurves += objcurves_content.from(objdata_pos+12);
  if (!write_file_from_String("objcurves_specific.py", specific_objcurves)) {
    warning("Unable to write objcurves_specific.py.");
    return false;
  }
  // 2. Call Blender with the specific Blender Python file.
  String blendercmd(blender_exec_path+" -b -P objcurves_specific.py");
  if (system(blendercmd)!=0) {
    warning("Blender call failed with command: "+blendercmd);
    return false;
  }
  String gzipcmd("gzip -9 "+objpath+".blend");
  if (system(gzipcmd)!=0) {
    warning("Gzip call failed with command:"+gzipcmd);
    return false;
  }
  // Delete temporary Blender Python file.
  // ...or don't...
  return true;
}

connectivity_graph::connectivity_graph(network & _n): net(&_n) {
  N = _n.PLLRoot<neuron>::length();
  X = new long[N];
  Y = new long[N];
  R = new long[N];
  n = new neuronptr[N];
  double maxradius = 0.0, r;
  PLL_LOOP_FORWARD(neuron,_n.PLLRoot<neuron>::head(),1) {
    r = e->Radius();
    if (r>maxradius) maxradius = r;
  }
  double circumpherence = 6.0*maxradius*((double) N);
  r = circumpherence / (2.0*M_PI);
  graph_R = FIGSCALED(r);
  double anglestep = (2.0*M_PI)/((double) N), angle = 0.0;
  int i = 0; double s;
  maxstrength = 0.0;
  PLL_LOOP_FORWARD(neuron,_n.PLLRoot<neuron>::head(),1) {
    n[i] = e;
    X[i] = FIGSCALED(r*cos(angle));
    Y[i] = FIGSCALED(r*sin(angle));
    R[i] = FIGSCALED(e->Radius());
    angle += anglestep;
    i++;
    PLL_LOOP_FORWARD_NESTED(connection,e->OutputConnections()->head(),1,c) {
      s = fabs(c->Abstract_Strength());
      if (s>maxstrength) maxstrength = s;
    }
  }
}

int connectivity_graph::find(neuron * m) {
  for (int i = 0; i<N; i++) if (n[i]==m) return i;
  return -1;
}
