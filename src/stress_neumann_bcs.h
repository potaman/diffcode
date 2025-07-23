/**
 * @file   stress_neumann_bcs.h
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Fri Mar 22 17:48:39 2013
 * 
 * @brief  This is a container to hold the boundary conditions for the stress solver.  * 
 * 
 */

#ifndef __stress_neumann_bcs_h__
#define __stress_neumann_bcs_h__

// libmesh includes 
#include "libmesh/id_types.h" 
#include "libmesh/function_base.h" 
#include "stress_solver.h" 

class DiffCode::StressSolver::StressNeumannBC{
 public: 

  /** 
   * This is the constructor for the stress neumann boundary condition
   * basically a holder for the boundary conditions
   * @param stress_solver 
   * @param boundary_id 
   * @param boundary_condition 
   * 
   * @return 
   */
  StressNeumannBC(StressSolver&  stress_solver,
		  boundary_id_type& boundary_id, 
		  FunctionBase<RealGradient>& boundary_condition); 


  /** 
   * Function to activate boundary condition 
   * 
   */

  void activate_boundary_condition(); 
  
  /** 
   * Function to deactivate boundary condition 
   * 
   */
  void deactivate_boundary_condition(); 

  /** 
   * 
   * Is bool active
   * 
   * @return 
   */
  bool is_active();


  RealGradient get_boundary_condition_value( const Point& p);
 private:
  /// The parent stress solver object. 
  StressSolver& _stress_solver; 
  
  /// the boundar id where the boundary condition is applied 
  boundary_id_type _boundary_id; 
  
  /// The actual value of the boundary condition.
  FunctionBase<RealGradient>& _boundary_condition; 

 
  bool _active; 


};




#endif 
