/**
 * @file   interfacial_average.C
 * @author Subramanya <ssadasiv@me-98-2.dhcp.ecn.purdue.edu>
 * @date   Sun Aug 25 18:00:22 2013
 * 
 * @brief  Averages to the nodes. in the interfacial regions. 
 *         I am going to try and do this diffusively. 
 *         The reason is that diffusive transformation 
 *         will get rid of the small bumps in the middle. 
 *         while keeping the value the same on the central line. 
 * 
 */


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
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/parameters.h"
#include "interfacial_average.h" 
#include "math.h"

using namespace libMesh; 

/// This is kind of weird. but this seems to make 
// sense. 
Real DiffCode::InterfacialAverage::_alpha;
Real DiffCode::InterfacialAverage::_beta; 
Real DiffCode::InterfacialAverage::_int_slope;

DiffCode::InterfacialAverage::InterfacialAverage(EquationSystems& eqSys, 
						   const std::string& name, 
						   const unsigned int number):
  LinearImplicitSystem(eqSys,name,number), 
  _diff_code(dynamic_cast<DiffCode&>(eqSys)) 
{
  std::string system_name = std::string("DiffCode::interfacial_average"); 
  std::string active_region_name= _diff_code.split_system_name(
						    system_name, 
						    name
						    ); 

  _active_region =(*_diff_code._active_regions.find(active_region_name)).second; 
  add_variable("var",FIRST,LAGRANGE,_active_region->get_subdomains()); 
  assemble_before_solve=false; 
}


void DiffCode::InterfacialAverage::compute_interfacial_average(ExplicitSystem& in_system, 
							  ExplicitSystem& out_system, 
							  const std::string& in_variable_name, 
							  const std::string& out_variable_name,
							  std::string type)  
{

  matrix->zero(); 
  rhs->zero(); 
  Number eps = _diff_code._eps; 
 
  std::string ch_system_name = _diff_code.form_system_name(
							   std::string("DiffCode::linear_cahn_hilliard"),
							   _active_region->get_name()
							   ); 
  
  
  TransientLinearImplicitSystem& ch = 
    _diff_code.get_system<TransientLinearImplicitSystem>(ch_system_name); 

  
  NumericVector<Number>& ch_vector =  *(ch.current_local_solution);
  NumericVector<Number>& in_vector =  *(in_system.current_local_solution); 
  

  const unsigned int var     = variable_number("var"); 
  const unsigned int in_var  = in_system.variable_number(in_variable_name); 
  const unsigned int phi_var = ch.variable_number("phi"); 


  const DofMap& dof_map = get_dof_map(); 
  const DofMap& dof_map_ch = ch.get_dof_map();  
  const DofMap& dof_map_in = in_system.get_dof_map(); 


  std::vector<unsigned int> dof_indices_phi; 
  std::vector<unsigned int> dof_indices_in; 
  std::vector<unsigned int> dof_indices; 


  Real small_number = _diff_code._small_number; 

  FEType fe_var = variable_type(var); 
  FEType fe_phi = ch.variable_type(phi_var); 
  FEType fe_in  = ch.variable_type(in_var); 

  AutoPtr<FEBase> fe_base_var(FEBase::build(_diff_code._dim, fe_var));
  AutoPtr<FEBase> fe_base_phi(FEBase::build(_diff_code._dim, fe_var));
  AutoPtr<FEBase> fe_base_in(FEBase::build(_diff_code._dim, fe_var));

  QGauss q_gauss(_diff_code._dim,fe_var.default_quadrature_order()); 

  fe_base_var->attach_quadrature_rule(&q_gauss);
  fe_base_phi->attach_quadrature_rule(&q_gauss); 
  fe_base_in->attach_quadrature_rule(&q_gauss); 

  const std::vector<Real>& JxW = fe_base_var->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_v = fe_base_var->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_v = fe_base_var->get_dphi(); 
  const std::vector<std::vector<RealGradient> >& dphi_ch = fe_base_phi->get_dphi(); 
  const std::vector<std::vector<Real> >& phi_ch = fe_base_phi->get_phi(); 
  const std::vector<std::vector<Real> >& phi_in = fe_base_in->get_phi(); 

  
  MeshBase::const_element_iterator el = 
    _diff_code._mesh->active_local_elements_begin(); 
  MeshBase::const_element_iterator el_end = 
    _diff_code._mesh->active_local_elements_end(); 

  DenseMatrix<Number> Ke; 
  DenseVector<Number> Fe; 
  
  //subdomain_id_type active_sub_id = _diff_code._mesh->get_id_by_name(_diff_code._active_set); 

  for(;el !=el_end; ++el)
    {
      const Elem* elem=*el; 
      subdomain_id_type sub_id = elem->subdomain_id(); 
      
      if(_active_region->subdomain_id_inclusion(sub_id))
	{
	  dof_map.dof_indices(elem,dof_indices); 
	  dof_map_ch.dof_indices(elem,dof_indices_phi,phi_var); 
	  dof_map_in.dof_indices(elem,dof_indices_in,in_var); 
	  
	  fe_base_var->reinit(elem); 
	  fe_base_phi->reinit(elem); 
	  fe_base_in->reinit(elem); 
	  
	  const unsigned int n_dofs_ch = dof_indices_phi.size(); 
	  const unsigned int n_dofs_in = dof_indices_in.size(); 
	  const unsigned int n_dofs   = dof_indices.size(); 
	  
	  Fe.resize(n_dofs); 
	  Ke.resize(n_dofs,n_dofs); 

	  for (unsigned int qp = 0; qp<q_gauss.n_points(); qp++) 
	    {
	      Number phi_qp = 0.0; 
	      Gradient grad_phi_qp(0.0,0.0,0.0); 
	      Number in_qp = 0.0; 
	     
	      for (unsigned int i=0; i<n_dofs_ch;i++) 
		{
		  grad_phi_qp+=dphi_ch[i][qp]*ch_vector(dof_indices_phi[i]); 
		  phi_qp += phi_ch[i][qp]*ch_vector(dof_indices_phi[i]); 
		}
	      for (unsigned int i=0; i<n_dofs_in;i++) 
		in_qp+=phi_in[i][qp]*in_vector(dof_indices_in[i]);

	      Number loc_qp = grad_phi_qp.size_sq(); 
	      Number heaviside_qp = (1+phi_qp)/2.0; 
	      Number loc_2_qp = (1-phi_qp*phi_qp);
	      

	      for (unsigned int i=0; i<n_dofs; i++) 
		for (unsigned int j=0; j<n_dofs;j++) 
		{
		  if (type == std::string("heaviside"))
		    Ke(i,j)+=JxW[qp]*phi_v[i][qp]*phi_v[j][qp];
		  else if (type == std::string("delta_1")) 
		    Ke(i,j)+=JxW[qp]*phi_v[i][qp]*phi_v[j][qp];
		  else if (type == std::string("delta")) 
		    Ke(i,j)+=JxW[qp]*phi_v[i][qp]*phi_v[j][qp];
		  else if (type == std::string("diffusive"))
		    Ke(i,j)+=JxW[qp]*(100*dphi_v[i][qp]*dphi_v[j][qp]
				      +(loc_qp/(eps*eps*eps))*phi_v[i][qp]*phi_v[j][qp]);
		}
	      for (unsigned int i=0; i<n_dofs;i++) 
		{
		  if (type == std::string("heaviside"))
		    Fe(i)+=JxW[qp]*phi_v[i][qp]*heaviside_qp*in_qp; 
		  else if (type == std::string("delta_1")) 
		    Fe(i)+=JxW[qp]*phi_v[i][qp]*loc_2_qp*in_qp; 
		  else if (type == std::string("delta")) 
		    Fe(i)+=JxW[qp]*phi_v[i][qp]*loc_qp*in_qp; 
		  else if (type == std::string("diffusive")) 
		    Fe(i)+=JxW[qp]*phi_v[i][qp]*(loc_qp/(eps*eps*eps))*in_qp;

		}
	    }
	  dof_map.constrain_element_matrix_and_vector(Ke,Fe,dof_indices); 
	  matrix->add_matrix(Ke,dof_indices); 
	  rhs->add_vector(Fe,dof_indices); 
	}
      

    }
  matrix->close(); 
  rhs->close();
  solve(); 
  
  const unsigned int sys_num = number(); 
  const unsigned int out_sys_num = out_system.number(); 
  
  const unsigned int var_num = variable_number("var"); 
  const unsigned int out_var_num = out_system.variable_number("mu"); 

  NumericVector<Number>& out_vector = *(out_system.current_local_solution); 
   
  MeshBase::const_node_iterator node_it = _diff_code._mesh->local_nodes_begin(); 
  const MeshBase::const_node_iterator node_it_end = _diff_code._mesh->local_nodes_end(); 
  
  for(; node_it!=node_it_end;node_it++) 
    {
      Node* node = *node_it; 
      unsigned int n_comp = node->n_comp(out_sys_num,out_var_num); 
      
      if(n_comp>0) 
	{
	  unsigned int dof_number_out = node->dof_number(out_sys_num,out_var_num,0); 
	  unsigned int dof_number = node->dof_number(sys_num,var_num,0);
	  //	  libMesh::out<<"here\t"<< solution->el(dof_number)<<std::endl;
	  out_vector.set(dof_number_out,solution->el(dof_number));//current_local_solution->el(dof_number)); 
	}
    }
  out_vector.close();
}

Real DiffCode::InterfacialAverage::compute_window(Real phi_value)
{
  Real eps = _diff_code._eps; 
  Real x = std::asin(phi_value); 
  Real alpha = _int_slope; 
  
  Real llimit = (PI/2.0)*(alpha - 1.0) ; 
  Real rlimit = (PI/2.0)*(1.0 - alpha) ; 

  
  if ( (x<= -PI/2.0) || (x >= PI/2.0) )
    return 0.0; 
  else if ((x>-PI/2.0) && (x < llimit))
    return 0.5*(1.0 - std::cos((PI+2*x)/alpha)); 
  else if ((x<PI/2.0) && (x > rlimit)) 
    return 0.5*(1.0 - std::cos((PI -2*x)/alpha)); 
  else 
    return 1.0; 
}

void DiffCode::InterfacialAverage::smooth_ch_chemical_potential()
{
  matrix->zero(); 
  rhs->zero(); 

  std::string ch_system_name = _diff_code.form_system_name(
							   std::string("DiffCode::linear_cahn_hilliard"),
							   _active_region->get_name()
							   ); 
  
  
  TransientLinearImplicitSystem& ch = 
    _diff_code.get_system<TransientLinearImplicitSystem>(ch_system_name); 

  NumericVector<Number>& ch_vector =  *(ch.current_local_solution);
  
  NumericVector<Number>& ch_sol_vector = *(ch.solution);

  // sol_copy=ch_vector;
  
     
  
  const unsigned int var = variable_number("var"); 
 
  const unsigned int phi_var = ch.variable_number("phi"); 
  const unsigned int mu_var = ch.variable_number("mu"); 

  const DofMap& dof_map = get_dof_map(); 
  const DofMap& dof_map_ch = ch.get_dof_map();  
  


  std::vector<unsigned int> dof_indices_phi; 
  std::vector<unsigned int> dof_indices_mu; 
  std::vector<unsigned int> dof_indices; 


  Real small_number = _diff_code._small_number; 

  FEType fe_var = variable_type(var); 
  FEType fe_phi = ch.variable_type(phi_var); 
  FEType fe_mu  = ch.variable_type(mu_var); 

  AutoPtr<FEBase> fe_base_var(FEBase::build(_diff_code._dim, fe_var));
  AutoPtr<FEBase> fe_base_phi(FEBase::build(_diff_code._dim, fe_var));
  AutoPtr<FEBase> fe_base_mu(FEBase::build(_diff_code._dim, fe_var));
  

  QGauss q_gauss(_diff_code._dim,fe_var.default_quadrature_order()); 

  fe_base_var->attach_quadrature_rule(&q_gauss);
  fe_base_phi->attach_quadrature_rule(&q_gauss); 
  fe_base_mu->attach_quadrature_rule(&q_gauss); 

  const std::vector<Real>& JxW = fe_base_var->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_v = fe_base_var->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_v = fe_base_var->get_dphi(); 
  const std::vector<std::vector<RealGradient> >& dphi_ch = fe_base_phi->get_dphi(); 
  const std::vector<std::vector<Real> >& phi_ch = fe_base_phi->get_phi(); 
  const std::vector<std::vector<Real> >& phi_mu = fe_base_mu->get_phi(); 

  
  MeshBase::const_element_iterator el = 
    _diff_code._mesh->active_local_elements_begin(); 
  MeshBase::const_element_iterator el_end = 
    _diff_code._mesh->active_local_elements_end(); 

  DenseMatrix<Number> Ke; 
  DenseVector<Number> Fe; 


  // subdomain_id_type active_sub_id = _diff_code._mesh->get_id_by_name(_diff_code._active_set); 

  for(;el !=el_end; ++el)
    {
      const Elem* elem=*el; 
      subdomain_id_type sub_id = elem->subdomain_id(); 

      if(_active_region->subdomain_id_inclusion(sub_id))
	{
	  dof_map.dof_indices(elem,dof_indices); 
	  dof_map_ch.dof_indices(elem,dof_indices_phi,phi_var); 
	  dof_map_ch.dof_indices(elem,dof_indices_mu,mu_var); 
	  
	  fe_base_var->reinit(elem); 
	  fe_base_phi->reinit(elem); 
	  fe_base_mu->reinit(elem); 
	  
	  const unsigned int n_dofs_ch = dof_indices_phi.size(); 
	  const unsigned int n_dofs_mu = dof_indices_mu.size(); 
	  const unsigned int n_dofs    = dof_indices.size(); 
	  
	  Fe.resize(n_dofs); 
	  Ke.resize(n_dofs,n_dofs); 
	  
	  for (unsigned int qp = 0; qp<q_gauss.n_points(); qp++) 
	    {
	      Number phi_qp  = 0.0; 
	      Number mu_qp = 0.0; 
	      for (unsigned int i=0; i<n_dofs_ch;i++) 
		{
		  phi_qp+=phi_ch[i][qp]*ch_vector(dof_indices_phi[i]);
		}
	      Real surf = compute_window(phi_qp) +small_number;
	      Real bulk = 1.0 - surf +small_number;



	      for (unsigned int i=0; i<n_dofs_mu;i++) 
		mu_qp+=phi_mu[i][qp]*ch_vector(dof_indices_mu[i]);

	      for (unsigned int i=0; i<n_dofs; i++) 
		for (unsigned int j=0; j<n_dofs;j++) 
		{
		  Ke(i,j)+=JxW[qp]*(bulk*_beta*dphi_v[i][qp]*dphi_v[j][qp]
				    +surf*_alpha*phi_v[i][qp]*phi_v[j][qp]);
		}
	      for (unsigned int i=0; i<n_dofs;i++) 
		{
		  
		  Fe(i)+=JxW[qp]*surf*_alpha*phi_v[i][qp]*mu_qp; 
		  
		}
	    }
	  dof_map.constrain_element_matrix_and_vector(Ke,Fe,dof_indices); 
	  matrix->add_matrix(Ke,dof_indices); 
	  rhs->add_vector(Fe,dof_indices); 
	}
      

    }
  matrix->close(); 
  rhs->close();
  solve(); 


  
  const unsigned int sys_num = number(); 
  const unsigned int out_sys_num = ch.number(); 
  
  
  MeshBase::const_node_iterator node_it = _diff_code._mesh->local_nodes_begin(); 
  const MeshBase::const_node_iterator node_it_end = _diff_code._mesh->local_nodes_end(); 
  
  for(; node_it!=node_it_end;node_it++) 
    {
      Node* node = *node_it; 
      unsigned int n_comp = node->n_comp(out_sys_num,mu_var); 
      
      if(n_comp>0) 
	{
	  unsigned int dof_number_out = node->dof_number(out_sys_num,mu_var,0); 
	  unsigned int dof_number = node->dof_number(sys_num,var,0);
	  
	  ch_sol_vector.set(dof_number_out,current_local_solution->el(dof_number));
	  
	}
    }
  //solution->close();
  ch_sol_vector.close();
  ch.update();



}




void DiffCode::InterfacialAverage::compute_interfacial_concentration()
{
  libMesh::out<<"computing interfacial concentration"<<std::endl;
 
  matrix->zero(); 
  rhs->zero(); 
  
  std::string ch_system_name = _diff_code.form_system_name(
							   std::string("DiffCode::linear_cahn_hilliard"),
							   _active_region->get_name()
							   ); 
  
  
  TransientLinearImplicitSystem& ch = 
    _diff_code.get_system<TransientLinearImplicitSystem>(ch_system_name); 
  
  std::string other_res_sys_name = _diff_code.form_system_name(
							       std::string("DiffCode::other_results"),
							       _active_region->get_name()
							       ); 
  
  ExplicitSystem& other_results = 
    _diff_code.get_system<ExplicitSystem>(other_res_sys_name); 
  
  NumericVector<Number>& ch_vector =  *(ch.current_local_solution);
  
  const unsigned int var = variable_number("var"); 
  const unsigned int phi_var = ch.variable_number("phi"); 
  const unsigned int mu_var = ch.variable_number("mu"); 
  
  const DofMap& dof_map = get_dof_map(); 
  const DofMap& dof_map_ch = ch.get_dof_map();  
  
  std::vector<unsigned int> dof_indices_phi; 
  std::vector<unsigned int> dof_indices_mu; 
  std::vector<unsigned int> dof_indices; 


  Real small_number = _diff_code._small_number; 

  FEType fe_var = variable_type(var); 
  FEType fe_phi = ch.variable_type(phi_var); 
  FEType fe_mu  = ch.variable_type(mu_var); 

  AutoPtr<FEBase> fe_base_var(FEBase::build(_diff_code._dim, fe_var));
  AutoPtr<FEBase> fe_base_phi(FEBase::build(_diff_code._dim, fe_var));
  AutoPtr<FEBase> fe_base_mu(FEBase::build(_diff_code._dim, fe_var));

  QGauss q_gauss(_diff_code._dim,fe_var.default_quadrature_order()); 

  fe_base_var->attach_quadrature_rule(&q_gauss);
  fe_base_phi->attach_quadrature_rule(&q_gauss); 
  fe_base_mu->attach_quadrature_rule(&q_gauss); 

  const std::vector<Real>& JxW = fe_base_var->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_v = fe_base_var->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_v = fe_base_var->get_dphi(); 
  const std::vector<std::vector<RealGradient> >& dphi_ch = fe_base_phi->get_dphi(); 
  const std::vector<std::vector<Real> >& phi_ch = fe_base_phi->get_phi(); 
  const std::vector<std::vector<Real> >& phi_mu = fe_base_mu->get_phi(); 
  
  MeshBase::const_element_iterator el = 
    _diff_code._mesh->active_local_elements_begin(); 
  MeshBase::const_element_iterator el_end = 
    _diff_code._mesh->active_local_elements_end(); 

  DenseMatrix<Number> Ke; 
  DenseVector<Number> Fe; 
  
  //subdomain_id_type active_sub_id = _diff_code._mesh->get_id_by_name(_diff_code._active_set); 

  for(;el !=el_end; ++el)
    {
      const Elem* elem=*el; 
      subdomain_id_type sub_id = elem->subdomain_id(); 
      
      if(_active_region->subdomain_id_inclusion(sub_id))
	{
	  dof_map.dof_indices(elem,dof_indices); 
	  dof_map_ch.dof_indices(elem,dof_indices_phi,phi_var); 
	  dof_map_ch.dof_indices(elem,dof_indices_mu,mu_var); 
	  
	  fe_base_var->reinit(elem); 
	  fe_base_phi->reinit(elem); 
	  fe_base_mu->reinit(elem); 
	  
	  const unsigned int n_dofs_ch = dof_indices_phi.size(); 
	  const unsigned int n_dofs_mu = dof_indices_mu.size(); 
	  const unsigned int n_dofs   = dof_indices.size(); 
	  
	  Fe.resize(n_dofs); 
	  Ke.resize(n_dofs,n_dofs); 
	  
	  for (unsigned int qp = 0; qp<q_gauss.n_points(); qp++) 
	    {
	      Number phi_qp  = 0.0; 
	      Number mu_qp = 0.0; 
	      for (unsigned int i=0; i<n_dofs_ch;i++) 
		{
		  phi_qp+=phi_ch[i][qp]*ch_vector(dof_indices_phi[i]);
		}
	      Real surf = compute_window(phi_qp) +small_number;
	      Real bulk = 1.0 - surf +small_number;
	      
	      for (unsigned int i=0; i<n_dofs_mu;i++) 
		mu_qp+=phi_mu[i][qp]*ch_vector(dof_indices_mu[i]);

	      for (unsigned int i=0; i<n_dofs; i++) 
		for (unsigned int j=0; j<n_dofs;j++) 
		{
		  Ke(i,j)+=JxW[qp]*(bulk*_beta*dphi_v[i][qp]*dphi_v[j][qp]
				    +surf*_alpha*phi_v[i][qp]*phi_v[j][qp]);
		}
	      for (unsigned int i=0; i<n_dofs;i++) 
		{
		  Fe(i)+=JxW[qp]*surf*_alpha*phi_v[i][qp]*exp(-mu_qp/5); 
		}
	    }
	  dof_map.constrain_element_matrix_and_vector(Ke,Fe,dof_indices); 
	  matrix->add_matrix(Ke,dof_indices); 
	  rhs->add_vector(Fe,dof_indices); 
	}
    }
  matrix->close(); 
  rhs->close();
  solve_and_constrain(); 
  
  const unsigned int sys_num = number(); 
  const unsigned int out_sys_num = other_results.number(); 
  const unsigned int out_var = other_results.variable_number("conc_int"); 

  NumericVector<Number>& out_vector = *(other_results.solution);
  
  MeshBase::const_node_iterator node_it = _diff_code._mesh->local_nodes_begin(); 
  const MeshBase::const_node_iterator node_it_end = _diff_code._mesh->local_nodes_end(); 
  
 
  for(; node_it!=node_it_end;node_it++) 
    {
      Node* node = *node_it; 
      unsigned int n_comp = node->n_comp(out_sys_num,mu_var); 
      unsigned int n_comp2 = node->n_comp(sys_num,var);
      if((n_comp>0) && (n_comp2 > 0)) 
	{
	  unsigned int dof_number_out = node->dof_number(out_sys_num,out_var,0);
	  unsigned int dof_number = node->dof_number(sys_num,var,0);
	  out_vector.set(dof_number_out,current_local_solution->el(dof_number));
	}
    }
  out_vector.close(); 
  other_results.update();
  
}

void DiffCode::InterfacialAverage::set_alpha(Real alpha)
{
  _alpha = alpha;
}

void DiffCode::InterfacialAverage::set_beta(Real beta)
{
  _beta = beta;
}

void DiffCode::InterfacialAverage::set_int_slope(Real int_slope)
{
  _int_slope = int_slope; 
}

void DiffCode::InterfacialAverage::solve_and_constrain()
{
  solve(); 
  get_dof_map().enforce_constraints_exactly(*this,solution.get());
  update();
}
