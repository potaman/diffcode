/**
 * @file   gradient_system.h
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Mon Mar 11 15:08:07 2013
 * 
 * @brief  This object holds the gradient field information at the nodes for all the 
 * scalar driving forces in the system. For the stress systems, this information is stored in a StressSystem object. 
 * 
 * 
 */

#ifndef __gradient_system_h__
#define __gradient_system_h__
//libmesh includes 
#include "libmesh/explicit_system.h" 

#include "diff_code.h"
#include "average_to_nodes.h"

class DiffCode::GradientSystem: public ExplicitSystem
{
  /** 
   * Constructor for the gradient system. Used to store all the gradient systems 
   * from the individual systems. the stress computations are stored in the stresssystem
   * @param eqSys  The parent equations systems object 
   * @param name  the name of this system 
   * @param number  the number of this system 
   * 
   * @return 
   */
 public:
  GradientSystem(
                 EquationSystems& eqSys, 
		   const std::string& name, 
		   const unsigned int& number
		 );

  /** 
   * Adds another field to compute the gradient of  
   * 
   * @param sys A reference to the system to which the 
   * @param var_name name of the variable whose gradients are to 
   *        stored here. 
   * Adds a variable with the name sys_name_var_name which is a first order lagrangian vector field.  
   */
  void AddField(
		const System& sys,
		const std::string& var_name
		);
 

  
  /** 
   * Function that calls the average to nodes object and averages all the fields
   * one by one. The function also has all defined as it's default argument, and will 
   * average all the gradient fields in this object by default. 
   * 
   * 
   * 
   */ 
  void ConstructGradientFields(); 
  /** 
   * Adds the number of the average to node system. 
   * 
   * @param a_to_n 
   */
  void add_average_to_node_system(const unsigned int  a_to_n); 
  
  
 private: 
  /// Reference to the parent diff_code object. 
  DiffCode& _diff_code; 
  /// variables for the system and the variable names. 
  std::multimap<std::string,std::string> _variables ; 


/// Average to nodes 
//  AverageToNodes& _a_to_n; 



};
#endif 
