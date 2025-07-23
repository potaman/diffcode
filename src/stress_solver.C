/**
 * @file   stress_solver.C
 * @author Subramanya <ssadasiv@me-98-66.dhcp.ecn.purdue.edu>
 * @date   Fri Jul  5 10:38:25 2013
 * 
 * @brief  this is the second attempt at the stress solver. 
 * the previous one did not work for some reason
 * 
 * 
 */
#include "stress_dirichlet_bcs.h" 
#include "stress_neumann_bcs.h"
#include "stress_solver.h"
#include "collated_phase_field.h"


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
#include "libmesh/dirichlet_boundaries.h"
	   
#include "libmesh/fe.h" 
	  

#include "libmesh/petsc_macro.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/petsc_linear_solver.h" 
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/parallel.h"

#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_implicit_system.h" 


#include "libmesh/const_function.h"
#include <sstream>

using namespace libMesh;

DiffCode::StressSolver::StressSolver(EquationSystems& eqSys, 
		       const std::string& name, 
		       const unsigned int number):
  LinearImplicitSystem(eqSys,name, number),
  _diff_code(dynamic_cast<DiffCode&>(eqSys))
{
  const unsigned int dim=_diff_code._dim; 
  if (dim>=2) 
    {
      unsigned int var = add_variable("u",FIRST,LAGRANGE); 
     
      var=add_variable("v",FIRST,LAGRANGE);
     
    }
  if (dim==3)
    {
     unsigned int var=add_variable("w",FIRST,LAGRANGE);
     
    }
  _constrained=false; 
}


void DiffCode::StressSolver::add_dirichlet_boundary_condition(boundary_id_type& i,
							      const std::vector<std::string> components, 
							      RealGradient value)
{

  
  std::set<boundary_id_type> bc_sides; 
  bc_sides.insert(i); 
  std::vector<unsigned int> variables(1);
  

  for (unsigned int i=0; i<components.size();i++)
    {
  
      
      if (components[i] == "x") 
	{
	  variables[0]=variable_number("u"); 
	  FunctionBase<Real>* bc_func = new ConstFunction<Real>(value(0)); 
	  DirichletBoundary* dbc = new DirichletBoundary(bc_sides, 
							 variables,
							 bc_func); 
	  get_dof_map().add_dirichlet_boundary(*dbc);
	  libMesh::out<<"constraining x"<<value(0)<<"\t";
	}
      else if (components[i] == "y") 
	{
	  variables[0]= variable_number("v"); 
	  FunctionBase<Real>* bc_func = new ConstFunction<Real>(value(1));
	  DirichletBoundary* dbc = new DirichletBoundary(bc_sides, 
							 variables,
							 bc_func); 
	  get_dof_map().add_dirichlet_boundary(*dbc);
	  libMesh::out<<"constraining y"<<value(1)<<std::endl;
	}
      else if ((components[i]=="z")&&(_diff_code._dim>2)) 
	{
	  
	  variables[0]= variable_number("w"); 
	  FunctionBase<Real>* bc_func = new ConstFunction<Real>(value(1));
	  DirichletBoundary* dbc = new DirichletBoundary(bc_sides, 
							 variables,
							 bc_func); 
	  get_dof_map().add_dirichlet_boundary(*dbc);
	  libMesh::out<<"constraining z";
	}
    }
}



void DiffCode::StressSolver::add_neumann_boundary_condition(boundary_id_type& i, 
							    RealGradient value)
{
  neumann_boundary_sides.insert(i); 
  FunctionBase<RealGradient>* bc_func =  new ConstFunction<RealGradient>(value); 

  StressNeumannBC* bc = new StressNeumannBC(*this, 
					    i, 
					    *bc_func); 

  _stress_neumann_bcs.insert(std::make_pair(i,*bc)); 
}





void DiffCode::StressSolver::reinit_lin_solver()
{
  // AutoPtr<PetscLinearSolver<Number> > petsc_solver = 
  //   AutoPtr<PetscLinearSolver<Number> > (new PetscLinearSolver<Number> ()); 

  PetscLinearSolver<Number>* petsc_solver = (PetscLinearSolver<Number>*) get_linear_solver(); 
 
  // Petsc error code object 
  PetscErrorCode ierr; 
  
  // Initialize the solver. creates the space  for the petsc
  petsc_solver->init();
  // Set the name of the petsc solver 

  KSP solver_ksp = petsc_solver->ksp();
  ierr =KSPSetOptionsPrefix(solver_ksp,"stress_solve_"); 
  
  // Set options from the options file 
  ierr =KSPSetFromOptions(solver_ksp); 
  
  // Put the solver back into the nonlinear solver 
  // linear_solver=petsc_solver; 
  release_linear_solver(petsc_solver);
}


void DiffCode::StressSolver::assemble() 
{
  const unsigned int dim = _diff_code._dim; 
  Number conc_0 =_diff_code.get_equilibrium_concentration(); 
  Number t_0 = _diff_code.get_equilibrium_temperature(); 
  
  
  Real small_number = _diff_code._small_number; 

  std::vector<unsigned int> var_numbers(dim); 
  
  if (dim >=2) 
    {
      var_numbers[0]=variable_number("u"); 
      var_numbers[1]=variable_number("v"); 
      
    }
  if (dim==3) 
    {
      var_numbers[2]=variable_number("w"); 
    }
  
 
  CollatedPhaseField& ch = *(_diff_code._collated_phase_field); 
 

  TransientLinearImplicitSystem& ts = _diff_code.get_system<TransientLinearImplicitSystem>("DiffCode::thermal");


  const NumericVector<Number>& ch_vector = *(ch.current_local_solution);
  const NumericVector<Number>& t_vector  = *(ts.current_local_solution);
  const unsigned int phi_var = ch.variable_number("phi"); 
  const unsigned int conc_var= ch.variable_number("conc");
  const unsigned int t_var = ts.variable_number("t"); 

  const DofMap& dof_map_str = get_dof_map(); 
  const DofMap& dof_map_ch = ch.get_dof_map(); 
  const DofMap& dof_map_t = ts.get_dof_map();
  
  
  FEType fe_ch  = ch.variable_type(phi_var); 
  FEType fe_str = variable_type(var_numbers[0]); // same variable numbers for everything. 
  FEType fe_t = ts.variable_type(t_var); 
  
  

  AutoPtr<FEBase> fe_base_str(FEBase::build(dim,fe_str)); 
  AutoPtr<FEBase> fe_base_str_surf(FEBase::build(dim,fe_str)); 
  
  AutoPtr<FEBase> fe_base_ch(FEBase::build(dim,fe_ch));
  AutoPtr<FEBase> fe_base_t(FEBase::build(dim,fe_t)); 


  QGauss q_bulk(dim,fe_str.default_quadrature_order()); 
  QGauss q_face(dim-1,fe_str.default_quadrature_order()); 

  
  fe_base_ch->attach_quadrature_rule(&q_bulk); 
  fe_base_str->attach_quadrature_rule(&q_bulk); 
  fe_base_str_surf->attach_quadrature_rule(&q_face);
  fe_base_t->attach_quadrature_rule(&q_bulk); 
  


  // Dof Indices for storing the degrees of freedom in.. 
  std::vector<unsigned int> dof_indices_str; 
  std::vector<std::vector<unsigned int> > dof_indices_str_comp(dim); 
  std::vector<unsigned int> dof_indices_t; 
 

  std::vector<unsigned int> dof_indices_ch_phi; 
  std::vector<unsigned int> dof_indices_ch_conc;
  // All the actual meat and potatoes of finite elements. 
  const std::vector<Real>& JxW = fe_base_str->get_JxW(); 
  const std::vector<Real>& JxW_surf = fe_base_str_surf->get_JxW(); 
 
  const std::vector<std::vector<Real> >& phi_str = fe_base_str->get_phi(); 
  const std::vector<std::vector<RealGradient> >& dphi_str = fe_base_str->get_dphi(); 
  const std::vector<std::vector<Real> >& phi_str_surf = fe_base_str_surf->get_phi(); 
  const std::vector<std::vector<RealGradient> >& dphi_str_surf = fe_base_str_surf->get_dphi(); 
  const std::vector<Point>& xyz_surf = fe_base_str_surf->get_xyz(); 
  const std::vector<Point>& normals_surf = fe_base_str_surf->get_normals(); 


  const std::vector<std::vector<Real> >& phi_ch = fe_base_ch->get_phi(); 
  const std::vector<std::vector<Real> >& phi_t = fe_base_t->get_phi(); 
  
  
  
  // Element iterator.. 
  MeshBase::const_element_iterator el =
    _diff_code._mesh->active_local_elements_begin(); 
  const MeshBase::const_element_iterator el_end = 
    _diff_code._mesh->active_local_elements_end(); 

  DenseVector<Number> Fe; 
  std::vector< DenseSubVector<Number>* > Fe_sub_vectors;
  //Fe_sub_vectors.resize(dim);
  for(unsigned int i=0; i< dim; i++) 
    {
      DenseSubVector<Number>* fe_sub = new DenseSubVector<Number>(Fe); 
      Fe_sub_vectors.push_back(fe_sub);
    }

  
  DenseMatrix<Number> Ke; 
  // The submatrices are arranged in rows. k11 k12 k13 k21 k22 k23 k31 k32 k33 .. 
  std::vector<DenseSubMatrix<Number>* > Ke_sub_matrices;
  //Ke_sub_matrices.resize();
  for (unsigned int i=0; i<dim*dim;i++) 
    {
      DenseSubMatrix<Number>* ke_sub = new DenseSubMatrix<Number>(Ke); 
      Ke_sub_matrices.push_back(ke_sub); 

    }

  DenseMatrix<Number> Km; 
  DenseMatrix<Number> Bl, Br; 
  //  set the appropriate sizes of the matrices. 

  if (dim==2) 
    {
      Km.resize(2,2); 
      Bl.resize(2,3); 
      Br.resize(2,3); 

    }
  else if (dim==3)
    {
      Km.resize(3,3); 
      Bl.resize(3,6); 
      Bl.resize(3,6); 
    }
    
  
  std::map<std::string,Material>::iterator mat_it;

  for(;el!=el_end;++el) 
    {
      const Elem* elem=*el;      
      subdomain_id_type sub_id = elem->subdomain_id();
      dof_map_str.dof_indices(elem,dof_indices_str); 
      
      for(unsigned int i=0; i<dim;i++) 
	{
	  dof_map_str.dof_indices(
				  elem,
				  dof_indices_str_comp[i],
				  var_numbers[i]
				  );
	}
      fe_base_str->reinit(elem); 
      fe_base_ch->reinit(elem);
      fe_base_t->reinit(elem);
      //get the degree of freedom for the thing.
      const unsigned int n_dofs = dof_indices_str.size();  
      const unsigned int n_dofs_comp = dof_indices_str_comp[0].size(); 
      
      
      // Resize to fit the current size of things. 
      Fe.resize(n_dofs); 
      Ke.resize(n_dofs,n_dofs); 
      
      // reposition the submatrices and the subvectors. 
      for (unsigned int i =0; i<dim; i++)
	{
	  for(unsigned int j=0; j < dim; j++) 
	    {
	      Ke_sub_matrices[dim*i+j]->reposition(var_numbers[i]*n_dofs_comp,
						   var_numbers[j]*n_dofs_comp,
						   n_dofs_comp,
						   n_dofs_comp); 
	    }
	  Fe_sub_vectors[i]->reposition(var_numbers[i]*n_dofs_comp,
				       n_dofs_comp); 
	}
      
      // Now check if this is an active subdomain and then use the 
      // information to compute the quantities
      std::string mat_name = ch.get_material_name(sub_id);

      
      mat_it = _diff_code._materials.find(mat_name); 
      
      
      DenseMatrix<Number> c_matrix = (mat_it->second).get_stiffness_matrix(); 
      unsigned int n_dofs_ch; 
      Number alpha = (mat_it->second).get_thermal_expansivity(); 
      Number omega = 0.0; //(mat_it->second).get_atomic_volume(); 
      

      if (ch.check_subdomain_inclusion(sub_id)) 
	{

	  dof_map_ch.dof_indices(elem,dof_indices_ch_phi,phi_var); 
	  dof_map_ch.dof_indices(elem,dof_indices_ch_conc,conc_var);
	  n_dofs_ch = dof_indices_ch_phi.size(); 
	  

	}else
	{
	  omega=0.0;
	}
      
      dof_map_t.dof_indices(elem,dof_indices_t); 
      unsigned int n_dofs_t =dof_indices_t.size(); 
      
      

      for (unsigned int qp=0; qp< q_bulk.n_points(); qp++) 
	{
	  Number heaviside = 1.0; 
	  Number conc = 0.0; 
	  if (ch.check_subdomain_inclusion(sub_id)) 
	    {
	      Real phi_qp_r = 0.0;
	      
	      for(unsigned int i=0; i< n_dofs_ch; i++) 
		{
		  conc += phi_ch[i][qp]*ch_vector(dof_indices_ch_conc[i]); 
		  phi_qp_r+=phi_ch[i][qp]*ch_vector(dof_indices_ch_phi[i]); 
		  
		}
	      heaviside = std::max((1+phi_qp_r)/2.0,0.0)+small_number; 
	     
	    }
	  Number t = 0.0; 
	  
	  for (unsigned int i = 0; i<n_dofs_t; i++) 
	    {
	      t+=phi_t[i][qp]*t_vector(dof_indices_t[i]); 
	    }
	 

	  

	  if (dim == 2) 
	    {
	      _update_2d_matrix(Ke_sub_matrices, 
				JxW,
				dphi_str,
				c_matrix,
				heaviside,
				n_dofs_comp,
				qp);
	      _update_2d_force_vector(Fe_sub_vectors,
				      JxW,
				      dphi_str,
				      c_matrix,
				      heaviside,
				      alpha*(t-t_0),
				      -omega*(conc-conc_0),
				      n_dofs_comp,
				      qp);

	    }
	  else 
	    {
	      _update_3d_matrix(Ke_sub_matrices, 
				JxW,
				dphi_str,
				c_matrix,
				heaviside,
				n_dofs_comp,
				qp);
	      _update_3d_force_vector(Fe_sub_vectors,
				      JxW,
				      dphi_str,
				      c_matrix,
				      heaviside,
				      alpha*(t-t_0),
				      -omega*(conc-conc_0),
				      n_dofs_comp,
				      qp);
	    }
	}
      
      for (unsigned int side=0; side<elem->n_sides(); side++) 
	if (elem->neighbor(side)==NULL) 
	  {
	  std::map<boundary_id_type,StressNeumannBC>::iterator it_nbc; 
	  
	  for(it_nbc=_stress_neumann_bcs.begin();
	      it_nbc!=_stress_neumann_bcs.end(); 
	      ++it_nbc)
	    {
	      if (_diff_code._mesh->boundary_info->has_boundary_id(elem,
								   side,
								   it_nbc->first))
		{

		  fe_base_str_surf->reinit(elem,side); 
		  StressNeumannBC* nbc = &(it_nbc->second); 
		  if (nbc->is_active()) 
		    {
		      for(unsigned int qp=0; qp<q_face.n_points();qp++) 
			{
			  RealGradient f = nbc->get_boundary_condition_value(xyz_surf[qp]);
			  for (unsigned int i=0; i<phi_str_surf.size(); i++) 
			    {
			      for(unsigned int j=0; j<dim;j++) 
				Fe_sub_vectors[j]->el(i)+=JxW_surf[qp]*(phi_str_surf[i][qp]*f(j));
			    }
			}
		    }
		}
	    }
	}
      dof_map_str.heterogenously_constrain_element_matrix_and_vector(Ke,Fe,dof_indices_str,false); 
      matrix->add_matrix(Ke,dof_indices_str); 
      rhs->add_vector(Fe,dof_indices_str);
    }
  matrix->close(); 
  rhs->close() ;
  
  //Call the destructor on all pointers in Ke_sub_matrics, Fe_sub_vectors;
  for(unsigned int i =0;i<dim*dim;i++) 
    {
      delete Ke_sub_matrices[i];
    }
  for(unsigned int i=0; i<dim; i++) 
    delete Fe_sub_vectors[i];

}


void DiffCode::StressSolver::_get_B_matrix(DenseMatrix<Number>& B, 
					   const  std::vector<std::vector<RealGradient> >& dphi,
					   const unsigned int i,
					   const unsigned int qp)
{
  if (_diff_code._dim == 2) 
    {
      B.resize(2,3);
      B(0,0) = dphi[i][qp](0); 
      B(1,1) = dphi[i][qp](1); 
      B(0,2) = dphi[i][qp](1);
      B(1,2) = dphi[i][qp](0);
    }
  else if (_diff_code._dim == 3) 
    {
      B.resize(3,6);
      B(0,0) = dphi[i][qp](0); 
      B(1,1) = dphi[i][qp](1);
      B(2,2) = dphi[i][qp](2);
      
      B(0,3) = dphi[i][qp](1);
      B(1,3) = dphi[i][qp](0);      

      B(1,4) = dphi[i][qp](2);
      B(2,4) = dphi[i][qp](1);
     
      B(0,5) = dphi[i][qp](2);
      B(2,5) = dphi[i][qp](0);

    }
}


void DiffCode::StressSolver::_get_stiffness_matrix_contrib(DenseMatrix<Number>& Km, 
							   DenseMatrix<Number>& Bl,
							   DenseMatrix<Number>& Br,
							   const DenseMatrix<Number>& stiffness, 
							   const Real& heaviside)
{
  Km*=0.0; 
  Km = Bl;
  Km.right_multiply(stiffness); 
  Km.right_multiply_transpose(Br);
  Km*=heaviside;
}




void DiffCode::StressSolver::_update_2d_matrix(std::vector<DenseSubMatrix<Number>* >& submatrices,
					       const std::vector<Real>& JxW,
					       const std::vector<std::vector<RealGradient> >& dphi,
					       const DenseMatrix<Number>& stiffness_matrix, 
					       const Real& heaviside,
					       const unsigned int& n_nodes, 
					       const unsigned int& qp)
{
  //unsigned int dim = 2;
  // the system submatrices are assumed to be arranged in the order. 
  // Kuu Kuv Kvu Kvv and so on . 
  // the stiffness matrix is a 3 by 3 matrix .. each of the Kuu , Kuv have been repositioned etc. 
  // I am pretty sure that part of the code is not slow.. it is the other stuff that is slowing the hell down. 
  for(unsigned int i=0; i< n_nodes; i++) 
    for(unsigned int j=0; j< n_nodes;j++) 
      {
	//Kuu
	submatrices[0]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](0)*stiffness_matrix(0,0)*dphi[j][qp](0)
						    +dphi[i][qp](1)*stiffness_matrix(2,2)*dphi[j][qp](1));
	//Kuv
	submatrices[1]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](0)*stiffness_matrix(0,1)*dphi[j][qp](1)
						    +dphi[i][qp](1)*stiffness_matrix(2,2)*dphi[j][qp](0));
	//Kvu
 	submatrices[2]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](1)*stiffness_matrix(1,0)*dphi[j][qp](0)
						    +dphi[i][qp](0)*stiffness_matrix(2,2)*dphi[j][qp](1)); 
	//Kvv
	submatrices[3]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](1)*stiffness_matrix(1,1)*dphi[j][qp](1)
						    +dphi[i][qp](0)*stiffness_matrix(2,2)*dphi[j][qp](0));
	
      }
}




void DiffCode::StressSolver::_update_3d_matrix(std::vector<DenseSubMatrix<Number>* >& submatrices,
					       const std::vector<Real>& JxW,
					       const std::vector<std::vector<RealGradient> >& dphi,
					       const DenseMatrix<Number>& stiffness_matrix, 
					       const Real& heaviside,
					       const unsigned int& n_nodes, 
					       const unsigned int& qp
					       )
{
  // The system submatrices are assumed to be arranged in ther order 
  // Kuu Kuv Kuw |  Kvu Kvv Kvw | Kwu Kwv Kww 



   for(unsigned int i=0; i< n_nodes; i++) 
    for(unsigned int j=0; j< n_nodes;j++) 
      {
	//Kuu

	submatrices[0]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](0)*stiffness_matrix(0,0)*dphi[j][qp](0)
						      +dphi[i][qp](1)*stiffness_matrix(5,5)*dphi[j][qp](1)
						      +dphi[i][qp](2)*stiffness_matrix(4,4)*dphi[j][qp](2)
						      );



	//Kuv 
	submatrices[1]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](0)*stiffness_matrix(0,1)*dphi[j][qp](1)
						      +dphi[i][qp](1)*stiffness_matrix(5,5)*dphi[j][qp](0));
						    
	
	
	
	// Kuw 
	submatrices[2]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](0)*stiffness_matrix(0,2)*dphi[j][qp](2)
						      +dphi[i][qp](2)*stiffness_matrix(4,4)*dphi[j][qp](0)); 
	

	
	// Kvu 
	submatrices[3]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](1)*stiffness_matrix(1,0)*dphi[j][qp](0)
						      +dphi[i][qp](0)*stiffness_matrix(5,5)*dphi[j][qp](1));
	
     
	// Kvv 
	submatrices[4]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](1)*stiffness_matrix(1,1)*dphi[j][qp](1)
						    +dphi[i][qp](2)*stiffness_matrix(3,3)*dphi[j][qp](2)
						    +dphi[i][qp](0)*stiffness_matrix(5,5)*dphi[j][qp](0));
	

	
	// Kvw
	submatrices[5]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](1)*stiffness_matrix(1,2)*dphi[j][qp](2)
						      +dphi[i][qp](2)*stiffness_matrix(3,3)*dphi[j][qp](1));

	
	// Kwu 
	submatrices[6]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](2)*stiffness_matrix(2,0)*dphi[j][qp](0)
						      +dphi[i][qp](0)*stiffness_matrix(4,4)*dphi[j][qp](2));

	
	
	// Kwv
	submatrices[7]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](2)*stiffness_matrix(2,1)*dphi[j][qp](1)
						      +dphi[i][qp](1)*stiffness_matrix(3,3)*dphi[j][qp](2));



	// Kww
	submatrices[8]->el(i,j) += heaviside*JxW[qp]*(dphi[i][qp](2)*stiffness_matrix(2,2)*dphi[j][qp](2)
						      +dphi[i][qp](0)*stiffness_matrix(4,4)*dphi[j][qp](0)
						      +dphi[i][qp](1)*stiffness_matrix(3,3)*dphi[j][qp](1));
	
      }







}


void DiffCode::StressSolver::_update_2d_force_vector(std::vector<DenseSubVector<Number>* >& subvectors,
						     const std::vector<Real>& JxW,
						     const std::vector<std::vector<RealGradient> >& dphi,
						     const DenseMatrix<Number>& stiffness_matrix,
						     const Real& heaviside,
						     const Real& deltaC,
						     const Real& deltaT,
						     const unsigned int& n_nodes,
						     const unsigned int& qp
						     )
{

  Number stress_1 = (stiffness_matrix(0,0) + stiffness_matrix(0,1))*(deltaT + deltaC);
  Number stress_2 = (stiffness_matrix(1,1) + stiffness_matrix(1,0))*(deltaT + deltaC); 
 
  for (unsigned int i=0; i < n_nodes; i++) 
    {
      //Fu
      subvectors[0]->el(i) += heaviside*JxW[qp]*(dphi[i][qp](0)*stress_1);
     
      //Fv
      subvectors[1]->el(i) += heaviside*JxW[qp]*(dphi[i][qp](1)*stress_2);
    }
}
			

void DiffCode::StressSolver::_update_3d_force_vector(std::vector<DenseSubVector<Number>* >& subvectors,
						     const std::vector<Real>& JxW,
						     const std::vector<std::vector<RealGradient> >& dphi,
						     const DenseMatrix<Number>& stiffness_matrix,
						     const Real& heaviside,
						     const Real& deltaC,
						     const Real& deltaT,
						     const unsigned int& n_nodes,
						     const unsigned int& qp
						     )
{
  Number stress_1 = (stiffness_matrix(0,0) + stiffness_matrix(0,1) + stiffness_matrix(0,2))*(deltaT+deltaC);
  Number stress_2 = (stiffness_matrix(1,0) + stiffness_matrix(1,1) + stiffness_matrix(1,2))*(deltaT+deltaC);
  Number stress_3 = (stiffness_matrix(2,0) + stiffness_matrix(2,1) + stiffness_matrix(2,2))*(deltaT+deltaC);



  for (unsigned int i=0; i< n_nodes; i++) 
    {
      //Fu
      subvectors[0]->el(i) += heaviside*JxW[qp]*(dphi[i][qp](0)*stress_1);
      //Fv
      subvectors[1]->el(i) += heaviside*JxW[qp]*(dphi[i][qp](1)*stress_2);
      //Fw
      subvectors[2]->el(i) += heaviside*JxW[qp]*(dphi[i][qp](2)*stress_3);
    }
}

void DiffCode::StressSolver::solve_and_constrain()
{
  solve(); 
  get_dof_map().enforce_constraints_exactly(*this,solution.get()); 
  update();
}


void DiffCode::StressSolver::free_matrix() 
{
  if (matrix->initialized())
    matrix->clear(); 
}

void DiffCode::StressSolver::allocate_matrix()
{
  if (!matrix->initialized()) 
    {
      matrix->init();
    }

}
