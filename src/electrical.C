/**
 * @file   electrical.C
 * @author Subramanya <subramanya@me-98-66.dhcp.ecn.purdue.edu>
 * @date   Tue Jun 11 10:53:29 2013
 * 
 * @brief  Implementation file for the electrical solver. Gets the phase field solver and then * and then uses the vector directly for everything.. This is a much simpler structure than 
 * whatever the hell I was trying earlier.. It is not possible to test this code without 
 * the electromigration part. 
 * 
 * 
 */
#include "libmesh/parameters.h"
#include "libmesh/boundary_info.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h" 
#include "libmesh/quadrature_gauss.h"
#include "libmesh/mesh.h"
#include "libmesh/fe.h" 
	  

#include "libmesh/petsc_macro.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/petsc_linear_solver.h" 
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/parallel.h"
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h" 
#include "libmesh/explicit_system.h"

#include "libmesh/const_function.h"


#include "electrical.h"
#include "collated_phase_field.h"
//=================================================================================
// Constructor for the electrical solver. 
//=================================================================================
DiffCode::Electrical::Electrical(EquationSystems& eqSys,
				 const std::string& name, 
				 const unsigned int number):
  LinearImplicitSystem(eqSys, name, number), 
  _diff_code(dynamic_cast<DiffCode&>(eqSys))
{
  const unsigned int dim = _diff_code._dim; 
  add_variable("e",FIRST,LAGRANGE); 
  _constrained=false; 
}

//=================================================================================
// End of constructor. 
//=================================================================================

//=================================================================================
// Adds a dirichlet boundary condition for the solver. 
//=================================================================================
void DiffCode::Electrical::add_dirichlet_boundary_condition(boundary_id_type& i, 
							       Real value)
{
  // Insert dirichlet boundary condition side into the solver. 
  _dirichlet_boundary_sides.insert(i); 
  
  // Create a function for the boundary condition... At present only constant. 
  FunctionBase<Real>* bc_func = new ConstFunction<Real>(value); 
  // Create the set of sides over which the current boundary condition needs to be applied. 
  std::set<boundary_id_type> bc_sides; 
  bc_sides.insert(i) ; 

  // Create the vector over which the boundary condition is to be applied. 
  std::vector<unsigned int> variables(1); 
  variables[0] = variable_number("e"); 

  // Create the boundary condition. 
  DirichletBoundary* dbc = new DirichletBoundary(bc_sides, 
						 variables, 
						 bc_func) ; 

  //Insert the boundary condition ito the set of boundary condiiton. 
  // this does not seem to be needed. need to figure out what to do . 
  _electrical_dirichlet_bcs.insert(std::make_pair(i,*dbc));
  get_dof_map().add_dirichlet_boundary((*dbc)); 
  // This sets that the code is constrained. 
  _constrained =true;
  libMesh::out<<"Added boundary condition"<<"\n";
  
}
//=================================================================================
// End of dirichlet boundary conditions. 
//=================================================================================

//=================================================================================
// Add the neumann boundary conditions. 
//=================================================================================
void DiffCode::Electrical::add_neumann_boundary_condition(boundary_id_type& i, 
							      RealGradient vector_value,
							      Real value ) 
  {
    // Insert the boundary condition sides into the neumann boundary condition set. 
    _neumann_boundary_sides.insert(i); 
    // Create the funcition for the boundary condition,
    FunctionBase<Real>* bc_func = new ConstFunction<Real>(value);
    // Create the function for the boundary condition
    FunctionBase<RealGradient>* bc_vec_func = new ConstFunction<RealGradient>(vector_value);   
    // Create a Neumann Boundary object
    NeumannBoundary* neumann_bc = new NeumannBoundary(i, 
						      *bc_vec_func,
						      *bc_func) ; 
    
    // Insert this into the map of boundary conditions. 
    _electrical_neumann_bcs.insert(std::make_pair(i,*neumann_bc));
    libMesh::out<<"Added boundary condition"<<"\n";
  }

//=================================================================================
// End of the neumann boundary condition
//=================================================================================

//=================================================================================
// Initializer for the linear solver. 
//=================================================================================

void DiffCode::Electrical::reinit_lin_solver()
{
  // AutoPtr<PetscLinearSolver<Number> > petsc_solver = 
  //   AutoPtr<PetscLinearSolver<Number> > (new PetscLinearSolver<Number> ()); 
 
  // Petsc error code object 
 
  PetscLinearSolver<Number>* petsc_solver =(PetscLinearSolver<Number>*) get_linear_solver(); 
  PetscErrorCode ierr; 
  
  // Initialize the solver. creates the space  for the petsc
  petsc_solver->init();
  // Set the name of the petsc solver 

  KSP solver_ksp = petsc_solver->ksp();
  ierr =KSPSetOptionsPrefix(solver_ksp,"elec_"); 
  
  // Set options from the options file 
  ierr =KSPSetFromOptions(solver_ksp); 
  
  // Put the solver back into the nonlinear solver 
  //linear_solver=petsc_solver; 

  release_linear_solver(petsc_solver); 
  
}


//=================================================================================
// Reinitializer for the linear solver. 
//=================================================================================

//=================================================================================
//Assembly routine for the assembly 
//=================================================================================
void DiffCode::Electrical::assemble()
{
  const unsigned int dim = _diff_code._dim; 
  
  Real small_number = _diff_code._small_number; 
  CollatedPhaseField& ch = *(_diff_code._collated_phase_field);
    
  const NumericVector<Number>& ch_vector = *(ch.current_local_solution); 

  
  const unsigned int e_var = variable_number("e"); 
  const unsigned int phi_var = ch.variable_number("phi"); 
  
  const DofMap& dof_map_e = get_dof_map(); 
  std::vector<unsigned int> dof_indices_e; 
    
  const DofMap& dof_map_ch = ch.get_dof_map(); 
  std::vector<unsigned int> dof_indices_ch_phi; 

  FEType fe_e = variable_type(e_var); 
  FEType fe_ch = ch.variable_type(phi_var); 

  AutoPtr<FEBase> fe_base_e(FEBase::build(dim,fe_e)); 
  AutoPtr<FEBase> fe_base_ch(FEBase::build(dim,fe_ch)); 

  QGauss q_bulk(dim, fe_e.default_quadrature_order()); 
  
  fe_base_e->attach_quadrature_rule(&q_bulk); 
  fe_base_ch->attach_quadrature_rule(&q_bulk); 

  const std::vector<Real>& JxW = fe_base_e->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_ch = fe_base_ch->get_phi();  
  const std::vector<std::vector<Real> >& phi_e = fe_base_e->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_e = fe_base_e->get_dphi(); 
  
  
  AutoPtr<FEBase> fe_base_surf_e(FEBase::build(dim,fe_e)); 
  QGauss q_face(dim-1,fe_e.default_quadrature_order()); 
 
  fe_base_surf_e->attach_quadrature_rule(&q_face); 
  
  const std::vector<Real>& JxW_surf = fe_base_surf_e->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_e_surf = fe_base_surf_e->get_phi(); 
  const std::vector<Point>& xyz_surf = fe_base_surf_e->get_xyz(); 
  const std::vector<Point>& norm_surf = fe_base_surf_e->get_normals(); 
  
  
  MeshBase::const_element_iterator el = 
    _diff_code._mesh->active_local_elements_begin(); 
  
  const MeshBase::const_element_iterator el_end = 
    _diff_code._mesh->active_local_elements_end(); 

  
  DenseVector<Number> Fe; 
  DenseMatrix<Number> Ke; 
  
 
  std::map<std::string,Material>::iterator mat_it;

  

  for(;el!=el_end;++el) 
    {
      // This solver needs to solve over the entire domain. 

      const Elem* elem = *el; 
      subdomain_id_type sub_id = elem->subdomain_id(); 
      
      dof_map_e.dof_indices(elem,dof_indices_e) ; 
      const unsigned int n_dofs = dof_indices_e.size(); 
      
      unsigned int n_dofs_ch; 
      

      
      if (ch.check_subdomain_inclusion(sub_id))
	{
	  dof_map_ch.dof_indices(elem, dof_indices_ch_phi, phi_var); 
	  n_dofs_ch = dof_indices_ch_phi.size(); 
	}
     
      Fe.resize(n_dofs); 
      Ke.resize(n_dofs,n_dofs); 

      fe_base_e->reinit(elem);  
      fe_base_ch->reinit(elem); 
      
      // compute the value of the phase field variables
      std::string mat_name = ch.get_material_name(sub_id); 

      mat_it = _diff_code._materials.find(mat_name); 
      Real conductivity = mat_it->second.get_conductivity();
     
      for (unsigned int qp=0; qp<q_bulk.n_points(); qp++)
 
	{
	  // If you are in the active subdomain. compute the phi_qp and then compute the conductivity
	  if(ch.check_subdomain_inclusion(sub_id))
	    {
	      Real phi_qp_r = 0.0; 
	      for (unsigned int i = 0;i< n_dofs_ch;i++) 
		phi_qp_r += phi_ch[i][qp]*ch_vector(dof_indices_ch_phi[i]); 
	      Real heaviside = std::max((1+phi_qp_r)/2.0,0.0) +small_number ; 
	      conductivity*=heaviside;
	    }
	  //Compute the matrix element
	  for (unsigned int i=0; i< n_dofs;i++) 
	      for(unsigned int j=0; j< n_dofs; j++) 
		Ke(i,j)+= JxW[qp]*conductivity*dphi_e[i][qp]*dphi_e[j][qp]; 
	  
	}
      // This is for the neumann boundary conditions. 
      for (unsigned int side =0; side<elem->n_sides(); side++) 
	{
	  if(elem->neighbor(side)==NULL)
	  {
	 
	    std::map<boundary_id_type,NeumannBoundary>::iterator it_nbc; 
	    
	    for (it_nbc=_electrical_neumann_bcs.begin(); 
		 it_nbc!=_electrical_neumann_bcs.end(); 
		 ++it_nbc) 
	      {
		
		//	libMesh::out<<it_nbc->first<<"here"; 
		//	libMesh::out<<it_nbc->second<<"\n";
		//      libMesh::out<<"\n";
		
		if(_diff_code._mesh->boundary_info->has_boundary_id(elem,side,it_nbc->first) )
		  {
		    fe_base_surf_e->reinit(elem,side); 
		    NeumannBoundary* nbc = &(it_nbc->second); 
		    if (nbc->is_active()) 
		      {
			fe_base_surf_e->reinit(elem,side); 
			for(unsigned int qp=0; qp< q_face.n_points(); qp++) 
			  {
			    Real val = nbc->get_boundary_condition_value(xyz_surf[qp], 
									 norm_surf[qp]); 
			    
			    for (unsigned int i =0; i < n_dofs ;i++)
			      Fe(i)+=JxW_surf[qp]*val*phi_e_surf[i][qp]; 
			  }
		      }
		  }
		
	      }
	  }
      	}
      dof_map_e.heterogenously_constrain_element_matrix_and_vector(Ke,Fe,dof_indices_e,false); 
      matrix->add_matrix(Ke, dof_indices_e); 
      rhs->add_vector(Fe, dof_indices_e); 
    }
  matrix->close();
  rhs->close(); 
}


void DiffCode::Electrical::solve_and_constrain()
{
  solve(); 
  get_dof_map().enforce_constraints_exactly(*this,solution.get());
  update();
  


}

void DiffCode::Electrical::compute_resistance(Number& current,
						Number& delta_v,
						Number& resistance,
						const std::string& side1,
						const std::string& side2) 
{
  current=0.0;
  delta_v=0.0;
  resistance=0.0;

  Number voltage_1=0.0;
  Number area_1=0.0; 
  Number voltage_2=0.0;
  Number area_2=0.0; 
  
  const unsigned int dim=_diff_code._dim; 
  const unsigned int e_var = variable_number("e"); 
  const DofMap& dof_map = get_dof_map(); 
  std::vector<unsigned int> dof_indices; 
    
  FEType fe = variable_type(e_var); 
  AutoPtr<FEBase> fe_base(FEBase::build(dim,fe)); 
  QGauss q_face(dim-1,fe.default_quadrature_order()); 
  fe_base->attach_quadrature_rule(&q_face); 
  
  const std::vector<Real>& JxW = fe_base->get_JxW(); 
  const std::vector<std::vector<Real> >& phi = fe_base->get_phi(); 
  const std::vector<std::vector<RealGradient> >& dphi = fe_base->get_dphi(); 
  const std::vector<Point>& norm = fe_base->get_normals(); 
  
  boundary_id_type s1 = _diff_code._mesh->boundary_info->get_id_by_name(side1);
  boundary_id_type s2 = _diff_code._mesh->boundary_info->get_id_by_name(side2); 

  
  MeshBase::const_element_iterator el = 
    _diff_code._mesh->active_local_elements_begin(); 
  const MeshBase::const_element_iterator el_end = 
    _diff_code._mesh->active_local_elements_end(); 
  
  
  for(;el!=el_end; ++el) 
    {
      const Elem* elem = *el; 
            
      dof_map.dof_indices(elem,dof_indices); 
      const unsigned int n_dofs = dof_indices.size(); 
    
      
      for (unsigned int side=0; side<elem->n_sides(); side++) 
	{
	  if (elem->neighbor(side)==NULL) 
	    {
	      bool s1t = _diff_code._mesh->boundary_info->has_boundary_id(elem,side,s1);
	      bool s2t = _diff_code._mesh->boundary_info->has_boundary_id(elem,side,s2);
	      
	      if (s1t) 
		{
		  fe_base->reinit(elem,side);
		  for (unsigned int qp=0; 
		       qp<q_face.n_points(); qp++) 
		    {
		      Gradient i_qp(0.0,0.0,0.0);
		      for (unsigned int i=0; i<n_dofs;i++)
			{
			  voltage_1+=JxW[qp]*phi[i][qp]*(current_local_solution->el(dof_indices[i]));
			  current+=JxW[qp]*(norm[qp]*dphi[i][qp])*(current_local_solution->el(dof_indices[i]));
			}
		    }
		}
	      else if (s2t)
		{
		  fe_base->reinit(elem,side); 
		  for (unsigned int qp=0; 
		       qp<q_face.n_points();qp++) 
		    {
		      for(unsigned int i=0; i<n_dofs; i++) 
			{
			  voltage_2+=JxW[qp]*phi[i][qp]*(current_local_solution->el(dof_indices[i]));
			}
		    }
		}
	    }
	}
    }
  
  libMesh::out<<"f"<<std::endl;
  _diff_code._mesh->comm().sum(voltage_1);
  _diff_code._mesh->comm().sum(voltage_2);
  _diff_code._mesh->comm().sum(current);
 
  delta_v = std::abs(voltage_1 - voltage_2);
  if(current > 1.e-14)
    {
      resistance = delta_v/std::abs(current);
    }
}
