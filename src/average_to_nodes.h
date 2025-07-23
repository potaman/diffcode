/**
 * @file   average_to_nodes.h
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Mon Mar 11 13:27:57 2013
 * 
 * @brief  Code to average gauss point fields to the nodes. 
 *     This is a function that allows to correctly and in parallel average gauss point fields to 
 *    the nodes without too much effort. 
 *    vector quantities like gradients and strains are averaged one component at a time. 
 *    A least squares matrix is computed to do the extrapolation. 
 * 
 */


#ifndef __average_to_nodes_h__
#define __average_to_nodes_h__ 


// libmesh includes 
#include "libmesh/linear_implicit_system.h" 
#include "diff_code.h" 
#include "material.h"
#include "stress_system.h"

class DiffCode::AverageToNodes: public LinearImplicitSystem, 
  public System::Assembly
{
 public: 
  /** 
   * Constructor for the average to nodes system 
   * 
   * @param eqSys the parent equation systems object.  
   * @param name  the name of this system. 
   * @param number the number for this system 
   * 
   * @return The Average to nodes. 
   */
  AverageToNodes(EquationSystems& eqSys, 
		 const std::string& name, 
		 const unsigned int number); 

  /** 
   * Assembly function for the object. 
   * 
   */
  void assemble(); 
  
  /// Boolean flag to indicate whether matrix is assembled 
  bool matrix_assembled; 

  /** 
   * Function to average the gradient of a vector variable to the  
   * nodes. Used primarily to compute the stress. 
   * @param out_system 
   * @param in_system 
   * @param in_variable_name 
   * @param out_variable_name
   */
  void average_gradient_to_nodes(			       
				 ExplicitSystem& out_system,
				 const ExplicitSystem& in_system, 
				 const std::string& in_variable_name,
				 const std::string& out_variable_name
				 );

  
 

  /** 
   * This is supposed to compute stresses and strains and put it back into the 
   * stress_system object. I tried to make it very general, but that doesn't 
   * work very well. 
   */ 
  void compute_stresses_and_strains( 
	
				     ); 


  
  /** 
   * This function smooths the chemical potential in the cahn hilliard problem. 
   * I am wondering whether to put it in the same place or to put it in a 
   * new results field.. 
   */
  void compute_interfacial_concentration(); 


  friend class StressSystem; 
  

  void solve_and_constrain(); 
  
  void set_diffusive_smoothing(Real alpha); 
  
  Real get_diffusive_smoothing(); 

 private:
  /// Pointer to the parent equation systems object. 
  DiffCode& _diff_code; 
  
  std::vector<std::string> _stress_variables; 
  std::vector<std::string> _strain_variables; 
  
  Real _diffusive_smoothing; 

  inline Real weird_heaviside(Real phi); 

};

#endif
