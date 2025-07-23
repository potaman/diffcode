/**
 * @file   gradient_system.C
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Tue Mar 12 01:01:49 2013
 * 
 * @brief  Implementation of gradient_system.h 
 * Stores the gradients of all the vectors as linear approximations of the vector at the nodes. 
 * This is then used to construct forcing vectors for the other systems. (Especially useful for the 
 * the gradients of the invariants. 
 */


#include "gradient_system.h"

#include <sstream>

DiffCode::GradientSystem::GradientSystem(
					 EquationSystems& eqSys, 
					 const std::string& name, 
					 const unsigned int& number					 ) : 
  ExplicitSystem(eqSys,name,number),
  _diff_code(dynamic_cast<DiffCode&>(eqSys))
{
  
}

void DiffCode::GradientSystem::AddField(
				   const System& sys, 
				   const std::string& var_name
				   )

{
  const std::string sys_name = sys.name(); 
  // add a pair to the multimap.. 
  _variables.insert(std::pair<std::string,std::string>(sys_name,
						       var_name)); 
  std::stringstream variable_name_stream ; 
  variable_name_stream<<sys_name<<"_"<<var_name;
  add_variable(variable_name_stream.str(),FIRST,LAGRANGE_VEC);
}


void DiffCode::GradientSystem::ConstructGradientFields( )
{
    
  std::multimap<std::string,std::string>::iterator it;

  for(it=_variables.begin(); it!=_variables.end();++it)
    {
      const ExplicitSystem& in_system = dynamic_cast<ExplicitSystem&>(_diff_code.get_system(it->first));
      std::stringstream outvar;
      outvar<<it->first<<"_"<<it->second;
      AverageToNodes& a_to_n = *(_diff_code._average_to_nodes);
      a_to_n.average_gradient_to_nodes( 
				       *this,
				       in_system, 
				       it->second, 
				       outvar.str()
					);
    }
}
