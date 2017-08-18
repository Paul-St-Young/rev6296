//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 */
template<class FT>
class DiffTwoBodyJastrowOrbital: public DiffOrbitalBase
{
  ///number of variables this object handles
  int NumVars;
  ///number of target particles
  int NumPtcls;
  ///number of groups, e.g., for the up/down electrons
  int NumGroups;
  ///variables handled by this orbital
  opt_variables_type myVars;
  ///container for the Jastrow functions  for all the pairs
  vector<FT*> F;
  ///offset for the optimizable variables
  vector<pair<int,int> > OffSet;
  Vector<RealType> dLogPsi;
  vector<GradVectorType*> gradLogPsi;
  vector<ValueVectorType*> lapLogPsi;
  std::map<std::string,FT*> J2Unique;
  ParticleSet* PtclRef;

public:

  ///constructor
  DiffTwoBodyJastrowOrbital(ParticleSet& p):NumVars(0)
  {
    PtclRef = &p;
    NumPtcls=p.getTotalNum();
    NumGroups=p.groups();
    F.resize(NumGroups*NumGroups,0);
  }

  ~DiffTwoBodyJastrowOrbital()
  {
    delete_iter(gradLogPsi.begin(),gradLogPsi.end());
    delete_iter(lapLogPsi.begin(),lapLogPsi.end());
  }

  void addFunc(int ia, int ib, FT* rf)
  { // add radial function rf to J2Unique
    SpeciesSet& species(PtclRef->getSpeciesSet());
    std::string pair_name = species.speciesName[ia] + species.speciesName[ib];
    J2Unique[pair_name]=rf;

    // also assign rf to correlate species pairs (ia,ib), (ib,ia)
    F[ia*NumGroups+ib]=rf;
    F[ib*NumGroups+ia]=rf; // enforce exchange symmetry
  }

  void linkFunc(int ia, int ib, FT* rf, bool exchange=true)
  { // assign FunctorType rf to species pair (ia,ib)
    F[ia*NumGroups+ib] = rf;
    if (exchange) F[ib*NumGroups+ia] = rf;
  }

  FT* getFunc(int ia, int ib)
  { // locate pair function by species indices (ia,ib)
    return F[ia*NumGroups+ib];
  }

  FT* findFunc(std::string pair_name)
  { // locate pair function by species pair name such as 'uu'
    FT* rf;
    bool found_coeff = false;
    typename std::map<std::string,FT*>::const_iterator
      it(J2Unique.begin()), it_end(J2Unique.end());
    for (;it!=it_end;it++)
    { // it->first: 'uu','ud'
      if (it->first==pair_name)
      {
        rf = it->second;
        found_coeff = true;
      }
    }
    if (!found_coeff) APP_ABORT(pair_name+" not found");
    return rf;
  }

  bool checkInitialization()
  { // make sure all pair functions stored in F are valid.
    //  namely, each functor pointer in F should exist in J2Unique, which is populated by addFunc
    bool all_found(true);
    for (int ifunc=0;ifunc<F.size();ifunc++)
    {
      FT* bsp = F[ifunc];
      // look for pair function in stored unique functions
      bool found(false);
      typename std::map<std::string,FT*>::const_iterator
        it(J2Unique.begin()),it_end(J2Unique.end());
      for (;it!=it_end;it++)
      {
        if (bsp == it->second) found=true;
      }
      if (!found) all_found = false;
    }
    return all_found;
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it++).second->resetParameters(active);
    }
  }

  ///reset the distance table
  void resetTargetParticleSet(ParticleSet& P)
  {
    PtclRef = &P;
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.clear();
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it).second->myVars.getIndex(active);
      myVars.insertFrom((*it).second->myVars);
      ++it;
    }
    myVars.getIndex(active);
    NumVars=myVars.size();

    //myVars.print(cout);

    if (NumVars && dLogPsi.size()==0)
    {
      dLogPsi.resize(NumVars);
      gradLogPsi.resize(NumVars,0);
      lapLogPsi.resize(NumVars,0);
      for (int i=0; i<NumVars; ++i)
      {
        gradLogPsi[i]=new GradVectorType(NumPtcls);
        lapLogPsi[i]=new ValueVectorType(NumPtcls);
      }
      OffSet.resize(F.size());
      int varoffset=myVars.Index[0];
      for (int i=0; i<F.size(); ++i)
      {
        if(F[i] && F[i]->myVars.Index.size())
        {
          OffSet[i].first=F[i]->myVars.Index.front()-varoffset;
          OffSet[i].second=F[i]->myVars.Index.size()+OffSet[i].first;
        }
        else
        {
          OffSet[i].first=OffSet[i].second=-1;
        }
      }
    }
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           vector<RealType>& dlogpsi,
                           vector<RealType>& dhpsioverpsi)
  {
    if(myVars.size()==0) return;
    bool recalculate(false);
    vector<bool> rcsingles(myVars.size(),false);
    for (int k=0; k<myVars.size(); ++k)
    {
      int kk=myVars.where(k);
      if (kk<0)
        continue;
      if (active.recompute(kk))
        recalculate=true;
      rcsingles[k]=true;
    }
    if (recalculate)
    {
      dLogPsi=0.0;
      for (int p=0; p<NumVars; ++p)
        (*gradLogPsi[p])=0.0;
      for (int p=0; p<NumVars; ++p)
        (*lapLogPsi[p])=0.0;
      vector<TinyVector<RealType,3> > derivs(NumVars);
      const DistanceTableData* d_table=P.DistTables[0];
      for (int i=0; i<d_table->size(SourceIndex); ++i)
      {
        for (int nn=d_table->M[i]; nn<d_table->M[i+1]; ++nn)
        {
          int ptype=d_table->PairID[nn];
          if(OffSet[ptype].first<0) continue; // nothing to optimize
          bool recalcFunc(false);
          for (int rcs=OffSet[ptype].first; rcs<OffSet[ptype].second; rcs++)
            if (rcsingles[rcs]==true) recalcFunc=true;
          if (recalcFunc)
          {
            std::fill(derivs.begin(),derivs.end(),0.0);
            if (!F[ptype]->evaluateDerivatives(d_table->r(nn),derivs)) continue;
            int j = d_table->J[nn];
            RealType rinv(d_table->rinv(nn));
            PosType dr(d_table->dr(nn));
            for (int p=OffSet[ptype].first, ip=0; p<OffSet[ptype].second; ++p,++ip)
            {
              RealType dudr(rinv*derivs[ip][1]);
              RealType lap(derivs[ip][2]+(OHMMS_DIM-1.0)*dudr);
              PosType gr(dudr*dr);
              dLogPsi[p]-=derivs[ip][0];
              (*gradLogPsi[p])[i] += gr;
              (*gradLogPsi[p])[j] -= gr;
              (*lapLogPsi[p])[i] -=lap;
              (*lapLogPsi[p])[j] -=lap;
            }
          }
        }
      }
      for (int k=0; k<myVars.size(); ++k)
      {
        int kk=myVars.where(k);
        if (kk<0)
          continue;
        if (rcsingles[k])
        {
          dlogpsi[kk]=dLogPsi[k];
          dhpsioverpsi[kk]=-0.5*Sum(*lapLogPsi[k])-Dot(P.G,*gradLogPsi[k]);
        }
        //optVars.setDeriv(p,dLogPsi[ip],-0.5*Sum(*lapLogPsi[ip])-Dot(P.G,*gradLogPsi[ip]));
      }
    }
  }

  DiffOrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    DiffTwoBodyJastrowOrbital<FT>* j2copy=new DiffTwoBodyJastrowOrbital<FT>(tqp);
    // 1. clone the map of unique pair functions "J2Unique"
    typename std::map<std::string,FT*>::const_iterator
      it(J2Unique.begin()), it_end(J2Unique.end());
    for (;it!=it_end;it++)
    {
      FT* rf = new FT(*it->second);
      j2copy->J2Unique[it->first] = rf;
    }

    // 2. map the array of pair functions "F"
    j2copy->F.resize(F.size());
    for (int i=0;i<F.size();i++)
    { // map unique functors into the non-unique functors,
      //  also check that every pair of particles have been assigned a functor
      bool done=false;
      typename std::map<std::string,FT*>::const_iterator it(j2copy->J2Unique.begin()),it_end(j2copy->J2Unique.end());
      for (;it!=it_end;it++)
      {
        done = true;
        j2copy->F[i] = it->second;
      }
      if (!done)
      {
        APP_ABORT("Error cloning DiffTwoBodyJastrowOrbital.\n")
      }
    }
    j2copy->myVars.clear();
    j2copy->myVars.insertFrom(myVars);
    j2copy->NumVars=NumVars;
    j2copy->NumPtcls=NumPtcls;
    j2copy->NumGroups=NumGroups;
    j2copy->dLogPsi.resize(NumVars);
    j2copy->gradLogPsi.resize(NumVars,0);
    j2copy->lapLogPsi.resize(NumVars,0);
    for (int i=0; i<NumVars; ++i)
    {
      j2copy->gradLogPsi[i]=new GradVectorType(NumPtcls);
      j2copy->lapLogPsi[i]=new ValueVectorType(NumPtcls);
    }
    j2copy->OffSet=OffSet;
    return j2copy;
  }

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1761 $   $Date: 2007-02-17 17:11:59 -0600 (Sat, 17 Feb 2007) $
 * $Id: TwoBodyJastrowOrbital.h 1761 2007-02-17 23:11:59Z jnkim $
 ***************************************************************************/

