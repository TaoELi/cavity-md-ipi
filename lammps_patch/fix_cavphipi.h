/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(cavphipi,FixCAVPHIPI)

#else

#ifndef LMP_FIX_CAVPHIPI_H
#define LMP_FIX_CAVPHIPI_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCAVPHIPI : public Fix {
 public:
  FixCAVPHIPI(class LAMMPS *, int, char **);
  virtual ~FixCAVPHIPI();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:
  char *host; int port; int inet, master, hasdata;
  int ipisock, me, socketflag; double *buffer; long bsize;
  // buffer for charge array
  double *buffer_q;
  int kspace_flag;
  int reset_flag;

 private:
  class Irregular *irregular;
};

}

#endif
#endif
