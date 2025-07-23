#ifndef __ch_neumann_boundary_h__
#define __ch_neumann_boundary_h__ 

#include "libmesh/system.h"
#include "libmesh/id_types.h" 

using namespace libMesh;

class CHNeumannBoundary
{
 public:
  /** 
   * Constructor
   * 
   * @param system  parent system where this boundary condition belongs to 
   * @param boundary_id  boundary where this boundary condition is to be applied
   * @param vector_condition  the vector part of this boundary condition
   * @param scalar_condition  the scalar part of this boundary condition 
   */
  CHNeumannBoundary(
		  const boundary_id_type& boundary_id, 
		  const Number& contact_angle,
		  const Number& surface_energy
		    ); 

  

  void activate_boundary_condition(); 

  void deactivate_boundary_condition(); 

  bool is_active(); 
  
  Number implicit_bc_rhs(const Number& phi_qp); 
  Number implicit_bc_lhs(const Number& phi_qp); 
  Number explicit_bc_rhs(const Number& phi_qp); 
 
  
  boundary_id_type get_boundary_id(); 
  
 private: 
  
  bool _active; 
  ///  Boundary where this boundary condition is applied
  boundary_id_type _boundary_id; 
  Number _surface_energy; 
  Number _contact_angle;
  
  
  
};
#endif // __neumann_boundary_h__




