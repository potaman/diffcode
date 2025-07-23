/**
 * @file   stress_dirichlet_bcs.C
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Fri Mar 22 18:01:44 2013
 * 
 * @brief  Implementation for the stress dirichlet boundary conditions. 
 * 
 * 
 */
#include "stress_dirichlet_bcs.h"



DiffCode::StressSolver::StressDirichletBC::StressDirichletBC(StressSolver& stress_solver, 
							     boundary_id_type& boundary_id, 
							     const std::vector<std::string>& constrained_coordinates, 
							     FunctionBase<RealGradient>& boundary_condition):
  _stress_solver(stress_solver), 
  _boundary_id(boundary_id), 
  _boundary_condition(boundary_condition)
{
  _constrained_coordinates.resize(3);
  for(unsigned int j=0; j<3;j++) 
    _constrained_coordinates[j]= false; 

  _constraint_vector*=0;
  for(unsigned int i=0; i<constrained_coordinates.size(); i++) 
    {

      if (constrained_coordinates[i]=="x") 
	{
	  _constrained_coordinates[0] = true;
	  _constraint_vector(0) = 1.0; 
	  
	}
      else if (constrained_coordinates[i]=="y") 
	{
	  _constrained_coordinates[1] = true; 
	  _constraint_vector(1) = 1.0;
	}
      else if (constrained_coordinates[i]=="z") 
	{
	  _constrained_coordinates[2] = true; 
	  _constraint_vector(2) = 1.0; 
	}

    }

  _active =true; 
}


void DiffCode::StressSolver::StressDirichletBC::activate_boundary_condition()
{
  _active = true; 
}

void DiffCode::StressSolver::StressDirichletBC::deactivate_boundary_condition()
{
  _active = false; 
}

void DiffCode::StressSolver::StressDirichletBC::unconstrain_coordinate(std::string unconstrained_coordinate)
{
  
  if (unconstrained_coordinate=="x") 
    {
      _constrained_coordinates[0] = false;
      _constraint_vector(0)*=0;
    }
  else if (unconstrained_coordinate=="y")
    {
      _constrained_coordinates[1] = false; 
      _constraint_vector(1)*=0;
    }
  else if (unconstrained_coordinate=="z") 
    {
      _constrained_coordinates[2] = false; 
      _constraint_vector(2)*=0;
    }

}


void DiffCode::StressSolver::StressDirichletBC::constrain_coordinate(std::string constrained_coordinate) 
{
 
  _constraint_vector*=0.0;
  if (constrained_coordinate=="x") 
    {
      _constrained_coordinates[0] = true;
      _constraint_vector(0) = 1.0; 
      
    }
  else if (constrained_coordinate=="y") 
    {
      _constrained_coordinates[1] = true; 
      _constraint_vector(1) = 1.0;
    }
  else if (constrained_coordinate=="z") 
    {
      _constrained_coordinates[2] = true; 
      _constraint_vector(2) = 1.0; 
    }
}

void DiffCode::StressSolver::StressDirichletBC::modify_boundary_condition(FunctionBase<RealGradient>& boundary_condition)
{

  _boundary_condition = boundary_condition; 
  
}

RealGradient DiffCode::StressSolver::StressDirichletBC::get_boundary_condition_value(const Point& xyz)
{
  return _boundary_condition(xyz);
}

RealGradient DiffCode::StressSolver::StressDirichletBC::get_constrained_coordinates()
{
  return _constraint_vector;
}



bool DiffCode::StressSolver::StressDirichletBC::is_active()
{
  return _active;
}



