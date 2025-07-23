/**
 * @file   initial_concentration.C
 * @author Subramanya <ssadasiv@me-98-66.dhcp.ecn.purdue.edu>
 * @date   Tue Jun 25 10:02:49 2013
 * 
 * @brief  Implementation of the diffusion solver for the system. 
 *         At present only considering the vacancy diffusion. 
 * 
 */


#include "initial_concentration.h" 




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
#include "libmesh/parallel.h"
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h" 
#include "libmesh/const_function.h"



//=================================================================================
// Constructor for the InitialConcentration  solver. 
//=================================================================================

DiffCode::InitialConcentration::InitialConcentration(EquationSystems& eqSys, 
			       const std::string& name, 
			       const unsigned int number):
  LinearImplicitSystem(eqSys,name,number), 
  _diff_code(dynamic_cast<DiffCode&>(eqSys))
{
  std::string system_name = std::string("DiffCode::init_conc"); 
  std::string active_region_name = _diff_code.split_system_name(
								system_name,
								name
						     ); 
  _active_region =(*_diff_code._active_regions.find(active_region_name)).second; 
  add_variable("init_conc",FIRST, LAGRANGE,_active_region->get_subdomains()); 
  attach_assemble_object(*this);
  _constrained=false; 
  _init_solve= true; 
}

///=================================================================================
/// Add dirichlet boundary condition for the solver. 
///=================================================================================
void DiffCode::InitialConcentration::add_dirichlet_boundary_condition(boundary_id_type& i, 
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
  variables[0] = variable_number("init_conc"); 

  // Create the boundary condition. 
  DirichletBoundary* dbc = new DirichletBoundary(bc_sides, 
						 variables, 
						 bc_func) ; 

  //Insert the boundary condition ito the set of boundary condiiton. 
  // this does not seem to be needed. need to figure out what to do . 
  _diffusion_dirichlet_bcs.insert(std::make_pair(i,*dbc));

  get_dof_map().add_dirichlet_boundary((*dbc)); 
  // This sets that the code is constrained. 
  _constrained =true;
  
}
//
//
//
//
///=================================================================================
/// Add neumann boundary condition 
///=================================================================================
//
void DiffCode::InitialConcentration::add_neumann_boundary_condition(boundary_id_type& i, 
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
   _diffusion_neumann_bcs.insert(std::make_pair(i,*neumann_bc));
  
 }


//
///=================================================================================
/// Assembly routine. 
///=================================================================================

void DiffCode::InitialConcentration::assemble() 
{
  rhs->zero();
  matrix->zero();
  const unsigned int dim = _diff_code._dim; 



  Real small_number = _diff_code._small_number; 
 


  std::string ch_sys_name = _diff_code.form_system_name("DiffCode::linear_cahn_hilliard",
							_active_region->get_name());

  
 
  
  TransientLinearImplicitSystem& ch = 
    _diff_code.get_system<TransientLinearImplicitSystem>(ch_sys_name); 
  

  const NumericVector<Number>& ch_vector = *(ch.current_local_solution); 
  const NumericVector<Number>& ch_vector_old = *(ch.old_local_solution); 
  
  const unsigned int conc_var = variable_number("init_conc"); 
  const unsigned int phi_var = ch.variable_number("phi"); 
  const unsigned int mu_var = ch.variable_number("mu"); 
  
  const DofMap& dof_map_conc = get_dof_map(); 
  const DofMap& dof_map_ch  = ch.get_dof_map(); 
 

  std::vector<unsigned int> dof_indices_conc; 
  std::vector<unsigned int> dof_indices_ch_phi; 
  std::vector<unsigned int> dof_indices_ch_mu; 
  


  FEType fe_conc = variable_type(conc_var); 
  FEType fe_ch = ch.variable_type(phi_var); 
  
  AutoPtr<FEBase> fe_base_conc(FEBase::build(dim, fe_conc)); 
  AutoPtr<FEBase> fe_base_ch(FEBase::build(dim,fe_ch)); 
  AutoPtr<FEBase> fe_base_surf_conc(FEBase::build(dim,fe_conc)); 
  
  QGauss q_bulk(dim,fe_conc.default_quadrature_order()); 
  QGauss q_face(dim-1,fe_conc.default_quadrature_order()); 

  
  fe_base_conc->attach_quadrature_rule(&q_bulk); 
  fe_base_ch->attach_quadrature_rule(&q_bulk); 
 
  fe_base_surf_conc->attach_quadrature_rule(&q_face); 



  const std::vector<Real>& JxW = fe_base_conc->get_JxW(); 
  const std::vector<std::vector<Real> >&  phi_conc = fe_base_conc->get_phi(); 
  const std::vector<std::vector<Real> >&  phi_ch = fe_base_ch->get_phi(); 
  

  const std::vector<std::vector<RealGradient> >& dphi_conc = fe_base_conc->get_dphi(); 
  const std::vector<std::vector<RealGradient> >& dphi_ch = fe_base_ch->get_dphi(); 
  

  const std::vector<Real>& JxW_surf = fe_base_surf_conc->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_conc_surf = fe_base_surf_conc->get_phi(); 
  const std::vector<Point>& xyz_surf = fe_base_surf_conc->get_xyz(); 
  const std::vector<Point>& norm_surf = fe_base_surf_conc->get_normals(); 


 
 
  MeshBase::const_element_iterator el = 
    _diff_code._mesh->active_local_elements_begin(); 
  const MeshBase::const_element_iterator el_end = 
    _diff_code._mesh->active_local_elements_end(); 
  
  std::map<boundary_id_type,NeumannBoundary>::iterator it_nbc; 
  DenseVector<Number> Fe; 
  DenseMatrix<Number> Ke; 
  
 
  // get an iterator for the material
  std::map<std::string,Material>::iterator mat_it; 
  
  // Find the material corresponding to the active set.
  // subdomain_id_type active_sub_id = _diff_code._mesh->get_id_by_name(_diff_code._active_set); 
  std::string mat_name = _active_region->get_material_name(); 
  
  mat_it = _diff_code._materials.find(mat_name);
  Material* mat = &(mat_it->second); 
  
  Real D = mat->get_diffusivity();  
  Real eps = _diff_code._eps; 
  Real dt  = _diff_code._dt;
         
  RealTensor diffT = mat->get_diffusivity_tensor(); 
  
  diffT.print();
 
  for (;el!= el_end;++el) 
    {
      const Elem* elem = *el; 
      subdomain_id_type sub_id = elem->subdomain_id(); 
     
      if (_active_region->subdomain_id_inclusion(sub_id)) 
	{
	  // Get the degrees of freedom.. 
	  dof_map_conc.dof_indices(elem,dof_indices_conc); 
	  dof_map_ch.dof_indices(elem,dof_indices_ch_phi,phi_var); 
	  dof_map_ch.dof_indices(elem,dof_indices_ch_mu,mu_var); 
	  

	  
	  fe_base_conc->reinit(elem); 
	  fe_base_ch->reinit(elem); 
	  
	  unsigned int n_dofs_conc   = dof_indices_conc.size(); 
	  unsigned int n_dofs_ch     = dof_indices_ch_phi.size(); 
	  

	  Ke.resize(n_dofs_conc,n_dofs_conc); 
	  Fe.resize(n_dofs_conc);
	  
	  for (unsigned int qp=0; qp<q_bulk.n_points(); qp++)
	    {
	      Real phi_qp = 0.0;
	      Real phi_old_qp = 0.0; 
	      RealGradient grad_phi_qp(0.0,0.0,0.0); 
	      RealGradient grad_phi_old_qp(0.0,0.0,0.0); 
	  
	      Real mu_qp = 0.0; 
	      Real mu_old_qp = 0.0; 
	      RealGradient grad_mu_qp(0.0,0.0,0.0); 
	      RealGradient grad_mu_old_qp(0.0,0.0,0.0); 


	      for (unsigned int i=0; i<n_dofs_ch; i++) 
		{
		  phi_qp+=phi_ch[i][qp]*ch_vector(dof_indices_ch_phi[i]);
		  phi_old_qp+=phi_ch[i][qp]*ch_vector_old(dof_indices_ch_phi[i]);

		  mu_qp+= phi_ch[i][qp]*ch_vector(dof_indices_ch_mu[i]); 
		  mu_old_qp+= phi_ch[i][qp]*ch_vector_old(dof_indices_ch_mu[i]); 

		  grad_phi_qp+=dphi_ch[i][qp]*ch_vector(dof_indices_ch_phi[i]); 
		  grad_phi_old_qp+=dphi_ch[i][qp]*ch_vector_old(dof_indices_ch_phi[i]); 
		  grad_mu_qp+=dphi_ch[i][qp]*ch_vector(dof_indices_ch_mu[i]); 
		  grad_mu_old_qp+=dphi_ch[i][qp]*ch_vector_old(dof_indices_ch_mu[i]); 
		}


	      Real delta = _delta_4_degree(phi_qp);
	      delta = (phi_qp>0.95)?0.0:delta;
	      Real heaviside = (1.0+phi_qp)/2.0 + small_number ;
	      //delta*=(std::abs(phi_qp)>0.9)?0.0:1.0; 
 		
	      for (unsigned int i=0;i<n_dofs_conc;i++) 
		for (unsigned int j=0; j<n_dofs_conc;j++)
		  {
		    Ke(i,j)+=JxW[qp]*(dphi_conc[i][qp]*dphi_conc[j][qp] 
						+((delta)/std::pow(eps,3.0))*phi_conc[i][qp]*phi_conc[j][qp]
				      );
		  }
	      for (unsigned int i=0; i<n_dofs_conc; i++) 
		{
		  
		  //Fe(i)+=JxW[qp]*((delta)/std::pow(eps,3.0))*phi_conc[i][qp]*std::exp(mu_qp);
		  Fe(i)+=JxW[qp]*((delta)/std::pow(eps,3.0))*phi_conc[i][qp];
		}
	    }
	  
	  // Go over the sides and apply the flux boundary conditions . I am hoping against hope
	  // that the thing can apply dirichlet boundary conditions over the internal 
	  // elements is possible. 
	  // For this one the neumann boundary condition is different. and instead of checking if the 
	  // adjacent element is null, we need to check if the adjacent element has the same subdomain name 
	  // as the current one. 
	  
	  for(unsigned int side=0; side<elem->n_sides(); side++) 
	    {
	      // Checking if the adjacent element does not belong to the same subdomain ID.. 
	      bool bc_check = (elem->neighbor(side)==NULL);
	      if (!bc_check)
		{
		  subdomain_id_type sub_id_neighbor = elem->neighbor(side)->subdomain_id()  ;
		 
		  bc_check = (sub_id != sub_id_neighbor) ; 
		}
	      if (bc_check) 
		{
		  std::map<boundary_id_type,NeumannBoundary>::iterator it_nbc; 
		  for(it_nbc=_diffusion_neumann_bcs.begin(); 
		      it_nbc!=_diffusion_neumann_bcs.end(); 
		      ++it_nbc)
		    {
		      if (_diff_code._mesh->boundary_info->has_boundary_id(elem,side,it_nbc->first))
			{
			  fe_base_surf_conc->reinit(elem,side);
			  NeumannBoundary* nbc = &(it_nbc->second); 
			  for(unsigned int qp = 0; qp<q_face.n_points(); qp++) 
			    {
			      Real val = nbc->get_boundary_condition_value(xyz_surf[qp], 
									   norm_surf[qp]); 
			      for (unsigned int i=0; i<n_dofs_conc; i++) 
				Fe(i)+=JxW_surf[qp]*val*phi_conc_surf[i][qp];
			      
			    }// for(unsigned int qp = 0; qp<q_face.n_points(); qp++) 

			}
		      
		    }// for(it_nbc=_diffusion_neumann_bcs.begin(); 
		  
		  
		}
	      
	      
	    }// for(unsigned int side=0; side<elem->n_sides(); side++) 
	  
	  dof_map_conc.heterogenously_constrain_element_matrix_and_vector(Ke,Fe,dof_indices_conc); 
	  matrix->add_matrix(Ke, dof_indices_conc); 
	  rhs->add_vector(Fe,dof_indices_conc); 
	  
	}
   }// for (;el!= el_end;++el) 
     
     
 matrix->close(); 
 rhs->close(); 

}
 
inline Number DiffCode::InitialConcentration::_delta_4_degree(const Number& phi)
{
  return std::max(std::pow((1.0 - phi)*(1.0 + phi),2),0.0);
}


