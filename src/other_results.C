/**
 * @file   other_results.C
 * @author Subramanya <ssadasiv@me-98-2.dhcp.ecn.purdue.edu>
 * @date   Tue Aug 13 11:28:34 2013
 * 
 * @brief  Implementation of other results.h
 * 
 * 
 */

#include "other_results.h" 


DiffCode::OtherResults::OtherResults(
				    EquationSystems& eqSys,
				    const std::string& name, 
				    const unsigned int  number
				    ):
  ExplicitSystem(eqSys,name,number),
  _diff_code(dynamic_cast<DiffCode&>(eqSys)) 
{
  std::string system_name = std::string("DiffCode::other_results"); 
  std::string active_region_name = _diff_code.split_system_name(
								system_name, 
								name
								); 

  _active_region=(*_diff_code._active_regions.find(active_region_name)).second; 
  
  add_variable("conc_int",FIRST,LAGRANGE,_active_region->get_subdomains());
}


