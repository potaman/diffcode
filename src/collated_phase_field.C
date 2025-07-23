/**
 * @file   collated_phase_field.C
 * @author Subramanya <ssadasiv@me-98-2.dhcp.ecn.purdue.edu>
 * @date   Thu Oct 17 13:58:39 2013
 * 
 * @brief  This is for the collated values of the 
 * phase field and the other things .
 * 
 * 
 */
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/error_vector.h"
#include "collated_phase_field.h" 
#include "active_region.h"

DiffCode::CollatedPhaseField::CollatedPhaseField(
						 EquationSystems& eqSys, 
						 const std::string& name, 
						 const unsigned int number
						 ):
  ExplicitSystem(eqSys,name,number),
  _diff_code(dynamic_cast<DiffCode&>(eqSys))
{
  std::map<std::string,ActiveRegion*>::iterator it;
  for (it=_diff_code._active_regions.begin(); it!=_diff_code._active_regions.end();++it) 
    {
      if(it->second->is_phase_field_region())
	_subdomain_ids.insert(it->second->get_subdomains()->begin(),
			      it->second->get_subdomains()->end());
      
      
      std::set<subdomain_id_type>::iterator sub_start = it->second->get_subdomains()->begin() ;
      const std::set<subdomain_id_type>::iterator sub_end = it->second->get_subdomains()->end(); 
      for(;sub_start!=sub_end; ++sub_start) 
	{

	  
	  _material_names.insert(std::pair<subdomain_id_type,std::string>(
									  *sub_start,
									  it->second->get_material_name()
									  )
				 );
}
    }
  add_variable("mu",FIRST,LAGRANGE,&_subdomain_ids); 
  add_variable("phi",FIRST,LAGRANGE,&_subdomain_ids); 
  add_variable("conc",FIRST,LAGRANGE,&_subdomain_ids); 
}
 
std::string DiffCode::CollatedPhaseField::get_material_name(subdomain_id_type sub_id) 
{
  
   if (_material_names.find(sub_id)!=_material_names.end()) 
     return _material_names.find(sub_id)->second; 
   else 
     return std::string("material property not found"); 
}
  
  bool DiffCode::CollatedPhaseField::check_subdomain_inclusion(subdomain_id_type sub_id)
{
return (_subdomain_ids.find(sub_id)!=_subdomain_ids.end());
}

    
void DiffCode::CollatedPhaseField::collate_phase_field_values()
{
/// A vector of the variable_numbers 
  std::vector<unsigned int> phi_var_numbers; 
  std::vector<unsigned int> mu_var_numbers; 
  std::vector<unsigned int> conc_var_numbers; 

  std::vector<TransientLinearImplicitSystem*> cahn_hilliard_systems; 
  std::vector<TransientLinearImplicitSystem*> diffusion_systems; 

  std::vector<unsigned int> diff_system_numbers; 
  std::vector<unsigned int> ch_system_numbers; 

  
  // // Go through the list of active regions 
  std::map<std::string,ActiveRegion*>::iterator it; 
  for (it=_diff_code._active_regions.begin(); it!=_diff_code._active_regions.end(); ++it)
    {
    
      if (it->second->is_phase_field_region()) 
   	{
   	  std::string diff_sys_name = _diff_code.form_system_name(
   								  std::string("DiffCode::diffusion"),
   								  it->first
   								  );
	  diffusion_systems.push_back(
				      &_diff_code.get_system<TransientLinearImplicitSystem>(diff_sys_name)
				      );


   	  std::string ch_sys_name = _diff_code.form_system_name(
  								std::string("DiffCode::linear_cahn_hilliard"),
   								it->first
   								);
	  cahn_hilliard_systems.push_back(
					  &_diff_code.get_system<TransientLinearImplicitSystem>(ch_sys_name) 
					  ); 
	  
   	}
    }

  for (unsigned int i=0; i< diffusion_systems.size(); i++) 
    {
      mu_var_numbers.push_back(cahn_hilliard_systems[i]->variable_number("mu"));
      phi_var_numbers.push_back(cahn_hilliard_systems[i]->variable_number("phi")); 
      ch_system_numbers.push_back(cahn_hilliard_systems[i]->number()); 
  

      conc_var_numbers.push_back(diffusion_systems[i]->variable_number("conc")); 
      diff_system_numbers.push_back(diffusion_systems[i]->number()); 
    }


  std::vector<NumericVector<Number>* > diff_vectors; 
  std::vector<NumericVector<Number>* > ch_vectors; 

  for (unsigned int i=0; i< diffusion_systems.size(); i++) 
    {
      diff_vectors.push_back((diffusion_systems[i]->current_local_solution.get())); 
      ch_vectors.push_back((cahn_hilliard_systems[i]->current_local_solution.get())); 
    }


   const unsigned int mu_var = variable_number("mu"); 
   const unsigned int phi_var = variable_number("phi"); 
   const unsigned int conc_var = variable_number("conc"); 
   
   MeshBase::const_node_iterator node_it = _diff_code._mesh->local_nodes_begin(); 
   const MeshBase::const_node_iterator node_it_end = _diff_code._mesh->local_nodes_end(); 

   unsigned int sys_num=number();  
   
   for(;node_it!=node_it_end;node_it++) 
     {
       Node* node = *node_it; 
       unsigned int n_dofs_mu = node->n_dofs(sys_num,mu_var); 
       unsigned int n_comp_mu = node->n_comp(sys_num,mu_var); 
       if(n_comp_mu > 0) 
	 {
	   for (unsigned int i=0; i <ch_system_numbers.size(); i++) 
	     if (node->n_comp(ch_system_numbers[i],mu_var_numbers[i]) > 0) 
	       {
		 unsigned int dof_number_mu = node->dof_number(sys_num,mu_var,0); 
		 unsigned int dof_number_mu_comp = node->dof_number(ch_system_numbers[i],mu_var_numbers[i],0); 
		 solution->set(dof_number_mu, 
			       ch_vectors[i]->el(dof_number_mu_comp)); 
		 unsigned int dof_number_phi = node->dof_number(sys_num, phi_var,0); 
		 unsigned int dof_number_phi_comp = node->dof_number(ch_system_numbers[i],phi_var_numbers[i],0); 
		 solution->set(dof_number_phi,
			       ch_vectors[i]->el(dof_number_phi_comp));
		 
		 unsigned int dof_number_conc = node->dof_number(sys_num,conc_var,0); 
		 unsigned int dof_number_conc_comp = node->dof_number(diff_system_numbers[i],conc_var_numbers[i],0); 
		 solution->set(dof_number_conc, 
			       diff_vectors[i]->el(dof_number_conc_comp)); 
	       }
	 }	
     }
   solution->close(); 
   update();
}


void DiffCode::CollatedPhaseField::refine_mesh_based_on_error()
{
  const unsigned int phi_var = variable_number("phi"); 
  const unsigned int mu_var = variable_number("mu"); 
  KellyErrorEstimator error_estimator; 
  ErrorVector error; 
  
  error_estimator.error_norm.set_weight(phi_var,1.0); 
  error_estimator.error_norm.set_weight(mu_var,0.0); 
  error_estimator.estimate_error(*this,error); 

  _diff_code._mesh_refinement->flag_elements_by_error_fraction(error);
 
  _diff_code._max_error_previous = _diff_code._max_error_current; 
  _diff_code._max_error_current  = error.maximum();  
  


}


void DiffCode::CollatedPhaseField::refine_mesh_based_on_values()
{

  const unsigned int phi_var = variable_number("phi"); 
  const unsigned int mu_var = variable_number("mu"); 
  const unsigned int max_levels = _diff_code._n_levels_refinement;
  // Store the current iterate in the vector  
  // const NumericVector<Number>& old_vector = *(old_local_solution); 
  const NumericVector<Number>& current_vector = *(current_local_solution); 


//   // Get the variables of the different systems
  FEType fe_type = variable_type(phi_var); 
  AutoPtr<FEBase> fe_phi (FEBase::build(_diff_code._dim,fe_type)); 
  QGauss qrule(_diff_code._dim,fe_type.default_quadrature_order());
  fe_phi->attach_quadrature_rule(&qrule); 
  
  const std::vector<std::vector<Real> >& phi_ch = fe_phi->get_phi(); 
  const DofMap& dof_map_ch = get_dof_map();
  
  std::vector<unsigned int> dof_indices_ch;
  std::vector<unsigned int> dof_indices_ch_phi; 
  
  MeshBase::const_element_iterator el = _diff_code._mesh->active_local_elements_begin(); 
  const MeshBase::const_element_iterator end_el = _diff_code._mesh->active_local_elements_end();  
  
  for(;el != end_el;++el)
    {
      Elem* elem = *el;
      // The assembly will only be done over the domain where the material 
      // is set as an active material 
      
      
      subdomain_id_type sub_id = elem->subdomain_id();
      if (check_subdomain_inclusion(sub_id)) 
	{
	  dof_map_ch.dof_indices(elem,dof_indices_ch); 
	  dof_map_ch.dof_indices(elem,dof_indices_ch_phi,phi_var); 
	  
	  const unsigned int n_dofs_ch_phi = dof_indices_ch_phi.size(); 

	  fe_phi->reinit(elem);
	  Number phi_av = 0.0; 
	  
	  for (unsigned int qp =0; qp<qrule.n_points();qp++)
	    {
	      for (unsigned int i=0; i<n_dofs_ch_phi;i++) 
		{
		  phi_av += phi_ch[i][qp]*current_vector(dof_indices_ch_phi[i]);
		}
	    }
	  phi_av/=qrule.n_points();
	  if ((phi_av>_diff_code._lower_pf_limit)&&(phi_av<_diff_code._upper_pf_limit)
	      && (elem->level() < max_levels))
	    {
	      elem->set_refinement_flag(Elem::REFINE);
	   
	    }
	  else
	    elem->set_refinement_flag(Elem::COARSEN);
	}
    }

}
