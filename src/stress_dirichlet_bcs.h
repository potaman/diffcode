/**
 * @file   stress_dirichlet_bcs.h
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Fri Mar 22 14:59:40 2013
 * 
 * @brief  This is just a container to hold the boundary conditions for the 
 *         stress solver. 
 * 
 * 
 */
#ifndef __stress_dirichlet_bcs_h__ 
#define __stress_dirichlet_bcs_h__ 


// libmesh includes 
#include "libmesh/id_types.h" 
#include "libmesh/function_base.h"
#include "stress_solver.h" 

class DiffCode::StressSolver::StressDirichletBC{
 public:
  
  /** 
   * Constructor for the StressDirichletBoundaryCondition  
   *  
   * @param stress_solver  stress solution object
   * @param boundary_id  boundary where it is to be constrained 
   * @param constrained_coordinates the coordinates that are to be constrained 
   * @param boundary_condition the function that applies the boundary condition. 
   * 
   * @return 
   */
  StressDirichletBC(StressSolver& stress_solver, 
		    boundary_id_type& boundary_id,
		    const std::vector<std::string>& constrained_coordinates, 
		    FunctionBase<RealGradient>& boundary_condition) ; 


  
  /** 
   * Constrain a single coordinate given the coordinate as a string 
   * 
   * @param constrained_coordinate 
   */
  void constrain_coordinate(std::string constrained_coordinate); 

  /** 
   * Unconstrain the coordinate given the coordinate as a string 
   * 
   * @param unconstrained_coordinate 
   */
  void unconstrain_coordinate(std::string unconstrained_coordinate); 
  

  /** 
   * Switches the boundary condition completely
   * 
   * @param boundary_condition 
   */
  void modify_boundary_condition(FunctionBase<RealGradient>& boundary_condition); 
 
  
  /** 
   * Activates the boundary condition 
   * 
   */
  void activate_boundary_condition(); 

  
  /** 
   * Deactivate the boundary condition 
   * 
   */

  void deactivate_boundary_condition(); 
 

  /** 
   * Returns the boundary condition value for the current 
   * 
   * @param xyz 
   * 
   * @return 
   */
  RealGradient get_boundary_condition_value(const Point& xyz); 
  
  
  /** 
   * Returns the active boundary condition at this point
   * 
   */
  RealGradient  get_constrained_coordinates();
    

  
  /** 
   * 
   * Returns whether the boundary condition is active. 
   * 
   * @return 
   */
  bool is_active()  ;    
    


  
 private:

  /// The parent stress solver object. 
  StressSolver& _stress_solver; 
  
  /// the boundary id where this boundary condition is applied
  boundary_id_type _boundary_id; 

  /// The boundary condition object. 
  FunctionBase<RealGradient>& _boundary_condition;
  

  /// Which of the coordinates are constrained 
  std::vector<bool> _constrained_coordinates; 

  
  RealGradient _constraint_vector;

  /// Check to see if the boundary condition in active. 
  bool _active; 
  
  

};




#endif 
