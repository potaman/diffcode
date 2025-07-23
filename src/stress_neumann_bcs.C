/**
 * @file   stress_neumann_bcs.C
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Fri Mar 22 20:00:45 2013
 * 
 * @brief  Neumann boundary conditions for the stress solver. implementation
 * 
 * 
 */

#include "stress_neumann_bcs.h"


DiffCode::StressSolver::StressNeumannBC::StressNeumannBC(StressSolver& stress_solver, 
							 boundary_id_type& boundary_id, 
							 FunctionBase<RealGradient>& boundary_condition):
  _stress_solver(stress_solver),
  _boundary_id(boundary_id),
  _boundary_condition(boundary_condition)
{
  _active=true; 
}


void DiffCode::StressSolver::StressNeumannBC::activate_boundary_condition()
{
  _active=true; 
}

void DiffCode::StressSolver::StressNeumannBC::deactivate_boundary_condition() 
{
  _active=false;
  
}

RealGradient DiffCode::StressSolver::StressNeumannBC::get_boundary_condition_value(const Point& p)
{
  return _boundary_condition( p);
}

bool DiffCode::StressSolver::StressNeumannBC::is_active()
{
  return _active;
}
