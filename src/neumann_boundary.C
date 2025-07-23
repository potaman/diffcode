/**
 * @file   neumann_boundary.C
 * @author Subramanya <subramanya@me-98-66.dhcp.ecn.purdue.edu>
 * @date   Mon Jun 10 13:42:47 2013
 * 
 * @brief  Implementation for neumann boundary conditions for 
 * most things that are not stress
 * 
 */


#include "neumann_boundary.h"


NeumannBoundary::NeumannBoundary( 
				 boundary_id_type& boundary_id, 
				 FunctionBase<RealGradient>& vector_condition, 
				 FunctionBase<Real>& scalar_condition):
  _boundary_id(boundary_id),
  _vector_condition(vector_condition),
  _scalar_condition(scalar_condition)
{
  _active=true; 

}


void NeumannBoundary::activate_boundary_condition()
{
  _active=true; 

}


void NeumannBoundary::deactivate_boundary_condition()
{
  _active=false; 


}

bool NeumannBoundary::is_active()
{
  return _active; 
}

Real NeumannBoundary::get_boundary_condition_value(const  Point& p , 
						   const Point& normal) 
{
  RealGradient temp; 
  Real value; 
  
  temp  = _vector_condition(p,0.0); 
  value = _scalar_condition(p,0.0); 
  
  for (unsigned int i=0; i<3;i++) 
    value+=(temp(i)*normal(i)); 
  
  return value;
}


boundary_id_type NeumannBoundary::get_boundary_id()
{
  return _boundary_id; 
  
}


