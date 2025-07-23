/**
 * @file   other_results.h
 * @author Subramanya <ssadasiv@me-98-2.dhcp.ecn.purdue.edu>
 * @date   Tue Aug 13 11:23:35 2013
 * 
 * @brief  This is an explicit system to hold all sorts of results for post processing. 
 * 
 * 
 */


#ifndef __other_results_h__
#define __other_results_h__ 

#include "libmesh/explicit_system.h" 
#include "diff_code.h" 
#include "average_to_nodes.h" 
#include "active_region.h"


class DiffCode::OtherResults: public ExplicitSystem
{
 public: 
  OtherResults(
	       EquationSystems& eqSys, 
	       const std::string& name, 
	       const unsigned int  number
	       ); 



  // void smooth_ch_chemical_potential(); 
  
 private: 
  DiffCode& _diff_code;
  ActiveRegion* _active_region; 
};
#endif
