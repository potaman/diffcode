
/**
 * @file   stress_solver.h
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Tue Mar 19 21:00:12 2013
 * 
 * @brief  Code to compute the stresses. 
 * 
 * 
 */
#ifndef __stress_solver_h__
#define __stress_solver_h__

#include "libmesh/linear_implicit_system.h" 
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/id_types.h"

#include "diff_code.h" 
#include "material.h"

class DiffCode::StressSolver:public LinearImplicitSystem,
  public System::Assembly
{
 public: 


  /** 
   * Constructor for the stress solver object. 
   * 
   * @param eqSys  parent equation systems object. 
   * @param name  the name of this system
   * @param number  the number of this system
   * 
   * @return 
   */
  StressSolver(EquationSystems& eqSys, 
	       const std::string& name, 
	       const unsigned int number) ; 

  /** 
   * Assembly function for the object. 
   * 
   */
  void assemble(); 

  /** 
   *  
   * Add a dirichlet boundary condition to the system
   * @param i 
   * @param components 
   * @param value 
   */
  void add_dirichlet_boundary_condition( boundary_id_type&  i,
					 const std::vector<std::string> components, 
					 RealGradient value
					 ); 
  
  
  /** 
   *  Add a neumann boundary condtion. 
   * 
   * @param i 
   * @param value 
   */
  void add_neumann_boundary_condition( boundary_id_type& i,
				       RealGradient value); 



  /** 
   * This is for reinitializing the solver. this I need to modify in the 
   * code 
   * 
   */
  void reinit_lin_solver();


  /// This frees memory
  void free_matrix(); 
  
  /// This allocates memory for the 
  /// matrix.. 
  void allocate_matrix(); 
  
  std::string _active_set; 
 

  void solve_and_constrain(); 

 private: 
  /// Pointer to the parent diff code object. 
  DiffCode& _diff_code;
  /// Reference to the materials object. 
 
  
  /// This is a class to contain the dirichlet boundary condition 
  friend class StressDirichletBC; 
  class StressDirichletBC; 

  /// This is a class to hold the dirichlet boundary condition 
  friend class StressNeumannBC; 
  class StressNeumannBC; 
  

  /// Listing of Dirichlet  boundary conditions 
  std::map<boundary_id_type,StressDirichletBC> _stress_dirichlet_bcs;
  /// Listing of Neumann boundary conditions 
  std::map<boundary_id_type,StressNeumannBC> _stress_neumann_bcs; 
  
  /// Dirichlet boundary sides. 
  std::set<boundary_id_type> dirichlet_boundary_sides; 
  
  /// Neumann Boundary Sides.
  std::set<boundary_id_type> neumann_boundary_sides; 
  
  /// To check if the boundary condition has been applied. 
  bool _constrained; 

  /// Materials 
  std::vector<Material> _mats;

  /** 
   * Compute the matrix contribution.  
   *
   * @param Km Densematrix to fill up the matrix  
   * @param Bl B matrix for the left side  
   * @param Br B matrix for the right side
   * @param sub subdomain id to find the material property 
   * @param param parameters for the evaluation of the material properties.
   */
  void _get_stiffness_matrix_contrib(DenseMatrix<Number>&  Km,
				     DenseMatrix<Number>&  Bl, 
				     DenseMatrix<Number>&  Br, 
				     const DenseMatrix<Number>& stiffness,
				     const Real&  heaviside );
  

  
  
  /** 
   * Get the B matrix. 
   * @param B 
   * @param dphi 
   * @param i ,
   * @param qp
   */
  void _get_B_matrix(DenseMatrix<Number>& B,
		     const std::vector<std::vector<RealGradient> >& dphi,
		     const unsigned int i,
		     const unsigned int qp); 
  



  /** 
   * Updates the sysstem submatrices for a 2d system
   * 
   * @param submatrices  Submatrices for the system 
   * @param JxW  Jacobian terms
   * @param dphi  Shape function derivatives in real space
   * @param stiffness_matrix  Stiffness matrix 
   * @param heaviside Heaviside function 
   * @param n_nodes number of nodes 
   * @param n_qp number of quadrature point
   */
  void _update_2d_matrix(std::vector<DenseSubMatrix<Number>* >& submatrices,
			 const std::vector<Real>& JxW,
			 const std::vector<std::vector<RealGradient> >& dphi,
			 const DenseMatrix<Number>& stiffness_matrix, 
			 const Real& heaviside,
			 const unsigned int& n_nodes, 
			 const unsigned int& n_qp
			 );
  
  
  /** 
   * Updates the system submatrices for a 3d system
   * 
   * @param submatrices 
   * @param JxW 
   * @param shape_functions 
   * @param stiffness_matrix 
   * @param Heaviside 
   * @param n_dofs 
   * @param qp 
   */
  void _update_3d_matrix(std::vector<DenseSubMatrix<Number>* >& submatrices,
			 const std::vector<Real>& JxW,
			 const std::vector<std::vector<RealGradient> >& dphi,
			 const DenseMatrix<Number>& stiffness_matrix, 
			 const Real& heaviside,
			 const unsigned int& n_nodes, 
			 const unsigned int& qp
			 );
  



  void _update_2d_force_vector(std::vector<DenseSubVector<Number>* >& subvectors,
			       const std::vector<Real>& JxW,
			       const std::vector<std::vector<RealGradient> >& dphi,
			       const DenseMatrix<Number>& stiffness_matrix,
			       const Real& heaviside,
			       const Real& deltaC,
			       const Real& deltaT,
			       const unsigned int& n_nodes,
			       const unsigned int& qp
			       );
			

  void _update_3d_force_vector(std::vector<DenseSubVector<Number>* >& subvectors,
			       const std::vector<Real>& JxW,
			       const std::vector<std::vector<RealGradient> >& dphi,
			       const DenseMatrix<Number>& stiffness_matrix,
			       const Real& heaviside,
			       const Real& deltaC,
			       const Real& deltaT,
			       const unsigned int& n_nodes,
			       const unsigned int& qp
			       );


  

 


};





#endif
