//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"


#include "Configuration.h"
#include "QMCHamiltonians/RPAPressure.h"
#include "QMCWaveFunctions/Jastrow/singleRPAJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "Optimize/VarList.h"
#include "ParticleBase/ParticleAttribOps.h"

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Message/Communicate.h"


namespace qmcplusplus
{

typedef QMCHamiltonianBase::Return_t Return_t;

void RPAPressure::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex=plist.add("Pressure");
  plist.add("ZVterm");
  plist.add("dpsi");
  plist.add("Edpsi");
  plist.add("Vdpsi");
}

void RPAPressure::setObservables(PropertySetType& plist)
{
  plist[myIndex]=Press;
  plist[myIndex+1]=Value;
  plist[myIndex+2]=tValue;
  plist[myIndex+3]=Energy*tValue;
  plist[myIndex+4]=Pot*tValue;
}

void RPAPressure::setParticlePropertyList(PropertySetType& plist, int offset)
{
  plist[myIndex+offset]=Press;
  plist[myIndex+1+offset]=Value;
  plist[myIndex+2+offset]=tValue;
  plist[myIndex+3+offset]=Energy*tValue;
  plist[myIndex+4+offset]=Pot*tValue;
}

void RPAPressure::resetTargetParticleSet(ParticleSet& P)
{
//     cout<<dPsi.size()<<endl;
  for (int i=0; i<(dPsi.size()); i++)
    dPsi[i]->resetTargetParticleSet(P);
  pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
  //3-D hard coded
  RealType tlen=std::pow(0.75/M_PI*P.Lattice.Volume/static_cast<RealType>(P.getTotalNum()),1.0/3.0);
  drsdV= tlen*pNorm;
//     app_log()<<"drsdV  "<<drsdV<<endl;
};

Return_t RPAPressure::evaluate(ParticleSet& P)
{
  vector<OrbitalBase*>::iterator dit(dPsi.begin()), dit_end(dPsi.end());
  tValue=0.0;
  Value=0.0;
  dG = 0.0;
  dL = 0.0;
  while(dit != dit_end)
  {
    tValue += (*dit)-> evaluateLog(P,dG,dL);
    ++dit;
  }
//     app_log()<<"tValue  "<<tValue<<endl;
  tValue *= drsdV;
  ZVCorrection =  -1.0 * (0.5*Sum(dL)+Dot(dG,P.G));
  Value = -ZVCorrection*drsdV;
  Pot=P.PropertyList[LOCALPOTENTIAL];
  Energy = P.PropertyList[LOCALENERGY];
  Press=2.0*Energy-Pot;
  Press*=pNorm;
  return 0.0;
}

Return_t RPAPressure::evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
{
  return evaluate(P);
}

bool RPAPressure::put(xmlNodePtr cur, ParticleSet& P)
{
  MyName = "RPAZVZBP";
  RealType tlen=std::pow(0.75/M_PI*P.Lattice.Volume/static_cast<RealType>(P.getTotalNum()),1.0/3.0);
  drsdV= tlen*pNorm;
  //Always a 2-body term
  //WARNING: NO CHIESA CORRECTION
  RPAJastrow* rpajastrow = new RPAJastrow(P,false);
  bool FOUND=false;
  xmlNodePtr tcur = cur->children;
  while(tcur != NULL)
  {
    string cname((const char*)tcur->name);
    if(cname == "TwoBody")
    {
      app_log() <<"  Found TwoBody parameters for dPsi in RPA-ZVZB Pressure"<<endl;
      rpajastrow->put(tcur);
      FOUND=true;
    }
    tcur = tcur->next;
  }
  if (!FOUND)
  {
    rpajastrow->buildOrbital("P_dRPA","true", "true", "drpa", tlen, 10.0/tlen);
    app_log() <<"  Using Default TwoBody parameters for dPsi in RPA-ZVZB Pressure"<<endl;
  }
  dPsi.push_back(rpajastrow);
  return true;
}

bool RPAPressure::put(xmlNodePtr cur, ParticleSet& P, ParticleSet& ISource, TrialWaveFunction& Psi)
{
  MyName = "HRPAZVZBP";
  RealType tlen=std::pow(0.75/M_PI*P.Lattice.Volume/static_cast<RealType>(P.getTotalNum()),1.0/3.0);
  drsdV= tlen*pNorm;
  //Always a 2-body term
  RPAJastrow* rpajastrow = new RPAJastrow(P,Psi.is_manager());
  bool FOUND1=false;
  bool FOUND2=false;
  xmlNodePtr tcur = cur->children;
  while(tcur != NULL)
  {
    string cname((const char*)tcur->name);
    if(cname == "TwoBody")
    {
      rpajastrow->put(tcur);
      dPsi.push_back(rpajastrow);
      app_log() <<"  Found TwoBody parameters for dPsi in RPA-ZVZB Pressure"<<endl;
      FOUND2=true;
    }
    else
      if(cname == "OneBody")
      {
        singleRPAJastrowBuilder* jb = new singleRPAJastrowBuilder(P, Psi, ISource);
        jb->put(tcur,0);
        OrbitalBase* Hrpajastrow = jb->getOrbital();
        dPsi.push_back(Hrpajastrow);
        app_log() <<"  Found OneBody parameters for dPsi in RPA-ZVZB Pressure"<<endl;
        FOUND1=true;
      }
    tcur = tcur->next;
  }
  if (!FOUND2)
  {
    rpajastrow->buildOrbital("P_dRPA","true", "true", "dRPA", tlen, 10.0/tlen);
    dPsi.push_back(rpajastrow);
    app_log() <<"  Using Default TwoBody parameters for dPsi in RPA-ZVZB Pressure"<<endl;
  }
//     if (!FOUND1) {
//       rpajastrow->buildOrbital("P_dRPA","true", "true", "dRPA", tlen, 10.0/tlen);
//       app_log() <<"  Using Default TwoBody parameters for dPsi in RPA-ZVZB Pressure"<<endl;
//     }
//     Communicate* cpt = new Communicate();
//     TrialWaveFunction* tempPsi = new TrialWaveFunction(cpt);
//     singleRPAJastrowBuilder* jb = new singleRPAJastrowBuilder(P, *tempPsi, source);
//     OrbitalBase* Hrpajastrow = jb->getOrbital();
//
//     delete jb;
//     delete tempPsi;
//     delete cpt;
//     dPsi.push_back(Hrpajastrow);
//     dPsi.push_back(rpajastrow);
  return true;
}

QMCHamiltonianBase* RPAPressure::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  RPAPressure* myClone = new RPAPressure(P);
  vector<OrbitalBase*>::iterator dit(dPsi.begin()), dit_end(dPsi.end());
  while(dit != dit_end)
  {
    myClone->dPsi.push_back((*dit)->makeClone(P));
    ++dit;
  }
  return myClone;
}

}


/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
* $Id: BareKineticEnergy.h 1581 2007-01-04 16:02:14Z jnkim $
***************************************************************************/

