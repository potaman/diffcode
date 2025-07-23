#ifndef __thermal_h__
#define __thermal_h__

#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/id_types.h"
#include "diff_code.h"
#include "material.h"
#include "neumann_boundary.h"




using libMesh::EquationSystems;
using libMesh::Gradient;
using libMesh::LinearImplicitSystem;
using libMesh::Number;
using libMesh::NumericVector;
using libMesh::Parameters;
using libMesh::Point;
using libMesh::SparseMatrix;
using libMesh::System;
using namespace::libMesh;


class DiffCode::Thermal: public TransientLinearImplicitSystem,
  public System::Assembly,
  public System::Initialization

{
 public:
  Thermal(EquationSystems& eqSys,
	  const std::string& name,
	  const unsigned int number);

  void assemble();

  void initialize(); 


  void solve_and_constrain();

  void reinit_lin_solver();

  void add_dirichlet_boundary_condition(boundary_id_type& i, 
					Real value); 


  void add_neumann_boundary_condition(boundary_id_type& i,
				      RealGradient vector_value,
				      Real value); 


  void add_convection_boundary_condition(boundary_id_type& i, 
					 Real co_eff,
					 Real ambient_temp); 

  std::string _active_set; 

  static Number initial_temp(const Point& p,
			     const Parameters& parameters,
			     const std::string&,
			     const std::string&);


 private:

  struct ConvectionBC{
    int bc_side; 
    Real ambient_temp; 
    Real co_eff;
  };

  bool _constrained; 
  std::map<boundary_id_type,DirichletBoundary> _thermal_dirichlet_bcs; 
  std::map<boundary_id_type,NeumannBoundary> _thermal_neumann_bcs; 
  std::map<boundary_id_type,ConvectionBC> _thermal_convection_bcs; 
  
  std::set<boundary_id_type> _dirichlet_boundary_sides; 
  std::set<boundary_id_type> _neumann_boundary_sides; 
  std::set<boundary_id_type> _convection_boundary_sides; 

  
  DiffCode& _diff_code; 
  
  
};



#endif
