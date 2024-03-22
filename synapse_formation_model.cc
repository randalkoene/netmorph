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
// synapse_formation_model.cc
// Randal A. Koene, 20050107

#include <math.h>
#include "diagnostic.hh"
#include "synapse_formation_model.hh"
#include "synapse_structure.hh"
#include "synapse.hh"
#include "network.hh"
#include "file.hh"

//#define DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS

synapse_formation_model defaultsfm(NULL,NULL);
synapse_formation_model * sfm = &defaultsfm;

general_synapse_formation_model_parameters_interface general_synapse_formation_model_parameters;

// An initial probability distribution of synapse types:
// Each value indicates the likelihood that a synaptic connection between
// the pre- and postsynaptic neurons includes a specific receptor type.
// A value less than 0 indicates that no such receptor can exist at that
// synaptic connection.
// [***INCOMPLETE] It would probably be better if arbor or terminal segments or
// spine objects had associated synapse formation models, which could determine
// locally what type of synaptic receptor channels and what type of synaptic
// connection can be established. This would be more flexible than a single
// multi-dimensional array.
double possible_syntypes[UNTYPED_NEURON+1][UNTYPED_NEURON+1][syntype_candidate] = {
  { // pre==PRINCIPAL_NEURON
    { // post==PRINCIPAL_NEURON
      1.0, -9.9, -9.9
    },
    { // post==INTERNEURON
      1.0, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      1.0, -9.9, -9.9
    },
    { // post==BIPOLAR
      1.0, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      1.0, 0.4, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==INTERNEURON
    { // post==PRINCIPAL_NEURON
      -9.9, -9.9, 1.0
    },
    { // post==INTERNEURON
      -9.9, -9.9, 1.0
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      -9.9, -9.9, 1.0
    },
    { // post==BIPOLAR
      -9.9, -9.9, 1.0
    },
    { // post==PYRAMIDAL
      -9.9, -9.9, 1.0
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==MULTIPOLAR_NONPYRAMIDAL
    { // post==PRINCIPAL_NEURON
      1.0, -9.9, -9.9
    },
    { // post==INTERNEURON
      1.0, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      1.0, -9.9, -9.9
    },
    { // post==BIPOLAR
      1.0, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      1.0, 0.4, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==BIPOLAR
    { // post==PRINCIPAL_NEURON
      1.0, -9.9, -9.9
    },
    { // post==INTERNEURON
      1.0, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      1.0, -9.9, -9.9
    },
    { // post==BIPOLAR
      1.0, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      1.0, 0.4, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==PYRAMIDAL
    { // post==PRINCIPAL_NEURON
      1.0, -9.9, -9.9
    },
    { // post==INTERNEURON
      1.0, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      1.0, -9.9, -9.9
    },
    { // post==BIPOLAR
      1.0, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      1.0, 0.9, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==UNTYPED_NEURON
    { // post==PRINCIPAL_NEURON
      -9.9, -9.9, -9.9
    },
    { // post==INTERNEURON
      -9.9, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      -9.9, -9.9, -9.9
    },
    { // post==BIPOLAR
      -9.9, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      -9.9, -9.9, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  }
};

// The following array provides pre- and postsynaptic neuron type specific
// information about the maximum distance between fibres at which a synapse
// can be established. Typically, a combination that commonly exhibits
// dendritic spikes can bridge a distance of up to 1.0 microns, while those
// without spines can bridge a distance of up to 0.1 microns.
double maxfibredistance[UNTYPED_NEURON+1][UNTYPED_NEURON+1] = {
  { // pre==PRINCIPAL_NEURON
      1.0, 0.1, 1.0, 1.0, 0.0
  },
  { // pre==INTERNEURON
      0.1, 0.1, 0.1, 0.1, 0.0
  },
  { // pre==MULTIPOLAR_NONPYRAMIDAL
      1.0, 0.1, 1.0, 1.0, 0.0
  },
  { // pre==PYRAMIDAL
      1.0, 0.1, 1.0, 1.0, 0.0
  },
  { // pre==UNTYPED_NEURON
      0.0, 0.0, 0.0, 0.0, 0.0
  }
};

uniform_pdf defaultsfmpdf(X_synapses);
probability_distribution_function * sfmpdf = &defaultsfmpdf;

bool no_autapses = true;
bool prev_conn_formed = false;


void general_synapse_formation_model_parameters_interface::parse_CLP(Command_Line_Parameters & clp) {
  // [***NOTE] This does not yet include a facility for setting these
  // parameters differently for specific synapse_formation_model objects.
  int n;
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1))
    for (neuron_type n_post = (neuron_type) 0; n_post <= UNTYPED_NEURON; n_post = (neuron_type) (n_post + 1))
      for (synapse_type receptortype = (synapse_type) 0; receptortype < syntype_candidate; receptortype = (synapse_type) (receptortype + 1))
	if ((n=clp.Specifies_Parameter("P_receptor."+String(neuron_short_name[n_pre])+'.'+String(neuron_short_name[n_post])+'.'+String(synapse_type_name[receptortype])))>=0) possible_syntypes[n_pre][n_post][receptortype] = atof(clp.ParValue(n));
  sfmpdf = pdfselection("synapse_formation",clp,X_synapses);
  if (!sfmpdf) sfmpdf = &defaultsfmpdf;
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1))
    for (neuron_type n_post = (neuron_type) 0; n_post <= UNTYPED_NEURON; n_post = (neuron_type) (n_post + 1))
      if ((n=clp.Specifies_Parameter("D_synmax."+String(neuron_short_name[n_pre])+'.'+String(neuron_short_name[n_post])))>=0) maxfibredistance[n_pre][n_post] = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("no_autapses"))>=0) no_autapses = (downcase(clp.ParValue(n))==String("true"));
}

String general_synapse_formation_model_parameters_interface::report_parameters() {
  String res("Synapse formation model hypotheses: ");
  res += "\n  Receptor type probabilities (P_receptor):\n                  POST\n";
  String setofreceptorsstr;
  for (synapse_type receptortype = (synapse_type) 0; receptortype < syntype_candidate; receptortype = (synapse_type) (receptortype + 1)) {
    String spacerstr("       ");
    for (int i = 0; ((i<6) && (synapse_type_name[receptortype][i]!='\0')); i++) spacerstr[i] = synapse_type_name[receptortype][i];
    setofreceptorsstr += spacerstr;
  }
  String setofneuronsstr;
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1)) {
#ifdef BIGSTRING_REPLICATE_PROBLEM
    int spacebufsize = (syntype_candidate*7)+2;
    char * spacebuf = new char[spacebufsize+1];
    for (int i = 0; i<spacebufsize; i++) spacebuf[i] = ' ';
    spacebuf[spacebufsize] = '\0';
    String spacerstr(spacebuf);
    delete[] spacebuf;
#else
    String spacerstr(replicate(' ',(syntype_candidate*7)+2));
#endif
    for (unsigned int i = 0; ((i<spacerstr.length()) && (neuron_short_name[n_pre][i]!='\0')); i++) spacerstr[i] = neuron_short_name[n_pre][i];
    setofneuronsstr += spacerstr;
  }
  res += "                  "+setofneuronsstr+'\n';
  res += "    PRE         ";
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1)) res += "  "+setofreceptorsstr;
  res += '\n';
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1)) {
    String spacerstr("                 ");
    for (int i = 0; ((i<12) && (neuron_short_name[n_pre][i]!='\0')); i++) spacerstr[i+4] = neuron_short_name[n_pre][i];
    res += spacerstr;
    for (neuron_type n_post = (neuron_type) 0; n_post <= UNTYPED_NEURON; n_post = (neuron_type) (n_post + 1)) {
      for (synapse_type receptortype = (synapse_type) 0; receptortype < syntype_candidate; receptortype = (synapse_type) (receptortype + 1)) res += String(possible_syntypes[n_pre][n_post][receptortype],"%6.2f ");
      if (n_post<UNTYPED_NEURON) res += "  ";
      else res += '\n';
    }
  }
  res += "  " + sfmpdf->report_parameters() + '\n';
  return res;
}

double synapse_formation_model::evaluate_possible_connection(fibre_segment * axonsegment, fibre_segment * dendritesegment) {
  // A model of synapse formation is applied to determine the likelihood that
  // a synapse forms between the axon segment and the dedrite segment.
  // In the base model, a non-zero likelihood leads to the creation of a
  // candidate synapse, i.e. a synapse with a development state set to the
  // virtual/quantum theoretical value of "candidate".
  // Returns the likelihood, which is also stored as a parameter of the
  // candidate synapse if one is created.
  // The calling function can decide to test connprocid in the other fibre
  // segment if it wishes to avoid evaluating the same pair of segments
  // more than once. Such multiple processing may otherwise occur in
  // specific cases, such as if segments that are allocated to multiple
  // spatial partitions appear in the same or neighboring paritions more
  // than once. In the current implementation of the growth model, this
  // does not otherwise occur, since each axon or dendrite segment is
  // evaluated only once for its possible connections with fibre segments
  // that developed earlier.
  // *** If multiple synapses between the same pair of segments are a
  // possibility then those can be created here or during a call to
  // establish_connections().
  DIAGNOSTIC_EVALUATE_POSSIBLE_CONNECTION_PROFILING;
  axonsegment->connprocid = dendritesegment;
  dendritesegment->connprocid = axonsegment;
  // Determine the likelihood
  if (no_autapses && (axonsegment->N()==dendritesegment->N())) return 0.0;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_synapse_structure);
  synapse_structure * s = new synapse_structure(*axonsegment,*dendritesegment); // ***WHY DO I CREATE THIS EVEN IF IT IS SUBSEQUENTLY NOT USED?
  DIAGNOSTIC_AFTER_ALLOCATION(new_synapse_structure,1,sizeof(*s));
#ifdef DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS
  double likelihood = 1.0 - (dist3D_Segment_to_Segment(*axonsegment,*dendritesegment,s)/proximitythreshold);
#else
double mfd = maxfibredistance[axonsegment->N()->TypeID()][dendritesegment->N()->TypeID()];
double likelihood = -1.0;


//"Crossing line-piece" synapse formation model:
//A geometrical approach to synapse formation. For a detailed explanation see
//SECO Year 2 report -- Grant agreement number 216593

double dsegseg = dist3D_Segment_to_Segment(*axonsegment,*dendritesegment,s);

//Initialise line-piece end points
double Px = 0.0;
double Py = 0.0;
double Pz = 0.0;
double Qx = 0.0;
double Qy = 0.0;
double Qz = 0.0;

double Rx = 0.0;
double Ry = 0.0;
double Rz = 0.0;
double Sx = 0.0;
double Sy = 0.0;
double Sz = 0.0;


double px = 0.0;
double py = 0.0;
double pz = 0.0;
double qx = 0.0;
double qy = 0.0;
double qz = 0.0;

double rx = 0.0;
double ry = 0.0;
double rz = 0.0;
double sx = 0.0;
double sy = 0.0;
double sz = 0.0;


fibre_segment * axon_test_seg;// = axonsegment;
fibre_segment * dendrite_test_seg;// = dendritesegment;

axon_test_seg = axonsegment;
dendrite_test_seg = dendritesegment;

#ifdef VECTOR3D
Px = axon_test_seg->P0.X();
Py = axon_test_seg->P0.Y();
Pz = axon_test_seg->P0.Z();
Qx = axon_test_seg->P1.X();
Qy = axon_test_seg->P1.Y();
Qz = axon_test_seg->P1.Z();


Rx = dendrite_test_seg->P0.X();
Ry = dendrite_test_seg->P0.Y();
Rz = dendrite_test_seg->P0.Z();
Sx = dendrite_test_seg->P1.X();
Sy = dendrite_test_seg->P1.Y();
Sz = dendrite_test_seg->P1.Z();
#endif



//I
//Perform spatial translation on axon and dendrite segments
//Relocate segments: Express their coordinates by shifting the
//x_0,y_0 point of the axon segment to the x=0,y=0,z=0 location.
px = Px - Px;
py = Py - Py;
pz = Pz - Pz;
qx = Qx - Px;
qy = Qy - Py;
qz = Qz - Pz;

rx = Rx - Px;
ry = Ry - Py;
rz = Rz - Pz;
sx = Sx - Px;
sy = Sy - Py;
sz = Sz - Pz;



//II
//First rotation: Rotate X2,Y2,Z2 coords of the axon segment into XY plane
//		  Rotate X1,Y1,Z1 and X2,Y2,Z2 coords of dendrite segment accordingly
double PQ = pow( (pow(abs(Px-Qx),2) + pow(abs(Py-Qy),2) + pow(abs(Pz-Qz),2) ), 0.5);
double l = pow( (pow(PQ,2) - pow((Qz - Pz),2)), 0.5);


double qx1;
double qy1;
double qz1;

qx1 = l;
qy1 = 0;
qz1 = qz;

double rx1;
double ry1;
double rz1;

rx1 = ((rx*qx) + (ry*qy))/l;
ry1 = ((-rx*qy) + (ry*qx))/l;
rz1 = rz;

double sx1;
double sy1;
double sz1;

sx1 = ((sx*qx) + (sy*qy))/l;
sy1 = ((-sx*qy) + (sy*qx))/l;
sz1 = sz;



//III
//Second rotation around Y-axis of q2 onto X-axis

double qx2;
double qy2;
double qz2;

qx2 = PQ;
qy2 = 0;
qz2 = 0;


double rx2;
double ry2;
double CR;
double rz2;

rx2 = ((rx*qx) + (ry*qy) + (rz*qz))/PQ;
ry2 = ry1;
CR = rx2;
rz2 = ((-CR*qz) + (PQ*rz))/l;

double sx2;
double CS;
double sy2;
double sz2;

sx2 = ((sx*qx)+(sy*qy)+(sz*qz))/PQ;
CS = sx2;
sy2 = ((-sx*qy)+(sy*qx))/l;
sz2 = ((-CS*qz) + (PQ*sz))/l;



//IV
//Third rotation around X-axis

double Or2_prime_sq;
double Os2_prime_sq;

double fv;

Or2_prime_sq = pow(ry2,2) + pow(rz2,2);
Os2_prime_sq = pow(sy2,2) + pow(sz2,2);
fv = 0.5*( Or2_prime_sq - Os2_prime_sq )/( pow((ry2 - sy2),2) + pow((rz2 - sz2),2) ) + 0.5;


if(abs(fv)<1e-10)
{
  fv = 0;
}

double Vy;
double Vz;
Vy = ry2 + fv*(sy2 - ry2);
Vz = rz2 + fv*(sz2 - rz2);

double lov;
lov = pow((pow(Vy,2) + pow(Vz,2)),0.5);

double rx3;
rx3 = CR;

double ry3;
ry3 = ((ry2*Vy) +(rz2*Vz))/lov;


double rz3;
rz3 = ((ry2*-Vz) +(rz2*Vy))/lov;

double sx3;
sx3 = CS;

double sy3;
sy3 = ((sy2*Vy) +(sz2*Vz))/lov;

double sz3;
sz3 = ((sy2*-Vz) +(sz2*Vy))/lov;



//V
//Determine whether r3s3 intersects XY-plane
bool CONDITION_A = false;
bool CONDITION_B1 = false;
bool CONDITION_B2 = false;
if(rz3*sz3 < 0)   //CONDITION A
{
  CONDITION_A = true;
}


double qx3 = PQ;


if(rx3 > sx3)   //CONDITION B1
{
  if(rx3 > 0 && sx3 < qx3)
  {
    CONDITION_B1 = true;
  }
}
else if(sx3 > rx3)
{
  if(sx3 > 0 && rx3 < qx3)
  {
    CONDITION_B1 = true;
  }
}

if(CONDITION_A == true && CONDITION_B1 ==true)
{
  double ux;  //Interception point
  ux = (CS*rz3 - CR*sz3)/(rz3-sz3);
  double tx;

  tx=ux;



  if(ux >=0 && ux<= PQ)  //Condidtion B2
  {
    CONDITION_B2 = true;
  }


  double frac_T_PQ = tx/qx;
  double frac_U_RS = (ux - rx3)/(sx3 - rx3);

  double Tx = Px + frac_T_PQ*(Qx - Px);
  double Ty = Py + frac_T_PQ*(Qy - Py);
  double Tz = Pz + frac_T_PQ*(Qz - Pz);

  double Ux = Rx + frac_U_RS*(Sx - Rx);
  double Uy = Ry + frac_U_RS*(Sy - Ry);
  double Uz = Rz + frac_U_RS*(Sz - Rz);

  double TU;
  TU= pow((pow((Tx - Ux),2) + pow((Ty - Uy),2) + pow((Tz - Uz),2)),0.5);
 
  
  if(mfd>=0.0)
  {
    //Set likelihood value for synapse formation between segments
    if(CONDITION_A == true && CONDITION_B1 == true && CONDITION_B2 == true && TU <= mfd)
    {
      likelihood = 1.0;
    }
    else
    {
      likelihood = 0.0;
    }
  }
}
else
{
  likelihood = 0.0;
}
#endif
  ///if ((eq->T()>=101000.0) && (eq->T()<102000.0)) {
  ///  cout << "\nSYNEVAL: T=" << eq->T() << " axonsegment=" << (unsigned long) axonsegment << " dendritesegment=" << (unsigned long) dendritesegment;
  ///}
  if (likelihood>0.0) {
    ///if ((eq->T()>=101000.0) && (eq->T()<102000.0)) cout << " CANDIDATE\n"; cout.flush();
    // Candidate synapses are created in the fibre structure, but not listed
    // in connections between neurons until the actual synapses are
    // established.
    connectionptr c = axonsegment->N()->connect_to(dendritesegment->N());
    new candidate_synapse(*c,likelihood,s);
    return likelihood;
    }
  ///if ((eq->T()>=101000.0) && (eq->T()<102000.0)) cout << '\n'; cout.flush();
  // plug a potential memory leak!
  DIAGNOSTIC_BEFORE_ALLOCATION(new_synapse_structure);
  delete s;
  DIAGNOSTIC_AFTER_ALLOCATION(new_synapse_structure,-1,sizeof(synapse_structure));
  return likelihood;
}

void synapse_formation_model::establish_connections() {
  // Selects which candidate synapses become actual synapses.
  // The type of actual neuron created depends on the location of the
  // synapse and on the types of neurons involved.
  if (!net) {
    warning("Warning: net==NULL in synapse_formation_model::establish_connections()\n");
    return;
  }
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    // ***synapse_inventory[syntype_GluR][1]++;
    PLL_LOOP_FORWARD_NESTED(connection,e->OutputConnections()->head(),1,c) {
      // ***synapse_inventory[syntype_iGluR][1]++;
      synapse * s_next;
      for (synapse * s = c->Synapses()->head(); (s); s = s_next) {
	s_next = s->Next(); // cached, since s may be removed
	// ***synapse_inventory[syntype_mGluR][1]++;
	if (s->type_ID()==syntype_candidate) {
	  // ***synapse_inventory[syntype_mGluR][0]++;
	  candidate_synapse * cs = static_cast<candidate_synapse *>(s);
	  // For this synapse formation model, the rules are outlined in
	  // Task Log entry TL#200508060747.1.
	  double likelihood = cs->Likelihood(); // == ratio of proximitythreshold microns
	  neuron_type pretype = cs->Presynaptic_Neuron()->TypeID();
	  neuron_type posttype = cs->Postsynaptic_Neuron()->TypeID();
#ifdef DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS
	  double distance = likelihood * proximitythreshold;
	  double mfd = maxfibredistance[pretype][posttype];
	  if (mfd<0.0) likelihood = -1.0;
	  else {
	    if (mfd==0.0) {
	      if (distance<0.00001) likelihood = 1.0;
	      else likelihood = 0.0;
	    }
	    else {
	      likelihood = 1.0 - (distance/mfd);
	    } 
	  }
	  if (likelihood>0.0) {
#endif
	    // ***synapse_inventory[syntype_GluR][0]++;
	    // Do the following for each possible synapse type (fully defined
	    // derived classes). This creates synapses (receptors) with initvalue
	    // ratio of initial amplitude parameter.
	    // [*** NOTE] The following method is based on an unsupported
	    // hypothesis with many assumptions that ignore the effects of
	    // early network activity. The method is described in
	    // TL#200508080920.5.
	    if (sfmpdf->random_positive()<=likelihood) { // a synapse is formed
	      candidates_conversion_attempted++;
	      double * initvalueptr = possible_syntypes[pretype][posttype];
	      if (X_synapses.get_rand_real1()<=initvalueptr[syntype_AMPAR]) new AMPA_receptors(*c,likelihood*INITIAL_STRENGTH_AMPA,cs->Structure());
	      if (X_synapses.get_rand_real1()<=initvalueptr[syntype_NMDAR]) new NMDA_receptors(*c,likelihood*INITIAL_STRENGTH_NMDA,cs->Structure());
	      if (X_synapses.get_rand_real1()<=initvalueptr[syntype_GABAR]) new GABA_receptors(*c,likelihood*INITIAL_STRENGTH_GABA,cs->Structure());
	    }
#ifdef DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS
	  }
#endif
	  // ***synapse_inventory[syntype_mGluR][0]++;
	  cs->remove(); // clean up candidate synapses
	}
      }
    }
  }
}
