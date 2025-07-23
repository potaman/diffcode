/**
 * @file   interfacial_average.h
 * @author Subramanya <ssadasiv@me-98-2.dhcp.ecn.purdue.edu>
 * @date   Sun Aug 25 17:04:59 2013
 * 
 * @brief  smooth interfacial quantities. 
 * 
 * 
 */


#ifndef __interfacial_average_h__
#define __interfacial_average_h__

#include "libmesh/linear_implicit_system.h"
#include "diff_code.h"
#include "active_region.h"


class DiffCode::InterfacialAverage: public LinearImplicitSystem
{
 public:
  InterfacialAverage(EquationSystems& eqSys, 
		     const std::string& name, 
		     const unsigned int number); 

  void compute_interfacial_average(ExplicitSystem& in_system,
				   ExplicitSystem& out_system, 
				   const std::string& in_variable_name, 
				   const std::string& out_variable_name,
				   std::string type); 

  

  void smooth_ch_chemical_potential(); 

  
  void compute_interfacial_concentration(); 

  /*
    set alpha for gradient smoothing
   */
  static void set_alpha(Real alpha);
  /*
    set beta for gradient smoothing
   */
  static void set_beta(Real beta);

  static void set_int_slope(Real int_slope); 

  void solve_and_constrain();
  

 private: 
  DiffCode& _diff_code; 
  static Real _alpha;
  static Real _beta;
  static Real _int_slope; 
  
  Real compute_window(Real phi_value); 
  
  ActiveRegion* _active_region;
  
  
};
#endif





