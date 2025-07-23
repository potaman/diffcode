/**
 * @file   stress_system.C
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Sat Mar 16 17:39:13 2013
 * 
 * @brief  Implementation of stress_system.h 
 * 
 */
#include "libmesh/explicit_system.h"
#include "libmesh/numeric_vector.h"
#include "stress_system.h"


#include <math.h>

using namespace libMesh; 

DiffCode::StressSystem::StressSystem(EquationSystems& eqSys, 
				     const std::string& name, 
				     const unsigned int& number
				     ):
  ExplicitSystem(eqSys,name,number),
  _diff_code(dynamic_cast<DiffCode&>(eqSys))
{
 
  // _diff_code= dynamic_cast<DiffCode&>(eqSys);
  // _a_to_n = a_to_n;
  
  const unsigned int dim = _diff_code._dim; 
  if(dim ==2) 
    {
      _stress_variable_numbers.resize(3); 
      _strain_variable_numbers.resize(3); 

      
      _stress_variable_numbers.push_back(add_variable("sigma_xx",FIRST,LAGRANGE));
      _stress_variable_numbers.push_back(add_variable("sigma_yy",FIRST,LAGRANGE));
      _stress_variable_numbers.push_back(add_variable("sigma_xy",FIRST,LAGRANGE));
      
      _stress_variables.resize(3); 
      _stress_variables[0]= "sigma_xx"; 
      _stress_variables[1]= "sigma_yy";
      _stress_variables[2] ="sigma_xy"; 

      
      _stress_variable_numbers.push_back(add_variable("epsilon_xx",FIRST,LAGRANGE)); 
      _stress_variable_numbers.push_back(add_variable("epsilon_yy",FIRST,LAGRANGE)); 
      _stress_variable_numbers.push_back(add_variable("epsilon_xy",FIRST,LAGRANGE)); 

      
      _strain_variables.resize(3); 
      _strain_variables[0]= "epsilon_xx"; 
      _strain_variables[1]= "epsilon_yy";
      _strain_variables[2]= "epsilon_xy"; 
 
    }

  else if (dim==3) 
    {
      _stress_variable_numbers.resize(6); 
      _strain_variable_numbers.resize(6); 
      
     _stress_variable_numbers.push_back( add_variable("sigma_xx",FIRST,LAGRANGE));
     _stress_variable_numbers.push_back( add_variable("sigma_yy",FIRST,LAGRANGE));
     _stress_variable_numbers.push_back( add_variable("sigma_zz",FIRST,LAGRANGE));
     _stress_variable_numbers.push_back( add_variable("sigma_xy",FIRST,LAGRANGE));
     _stress_variable_numbers.push_back( add_variable("sigma_xz",FIRST,LAGRANGE));
     _stress_variable_numbers.push_back( add_variable("sigma_yz",FIRST,LAGRANGE));
      
     _stress_variables.resize(6); 
     _stress_variables[0]= "sigma_xx"; 
     _stress_variables[1]= "sigma_yy";
     _stress_variables[2] ="sigma_zz"; 
     _stress_variables[3]= "sigma_xy"; 
     _stress_variables[4]= "sigma_xz";
     _stress_variables[5] ="sigma_yz"; 

      
      _strain_variable_numbers.push_back(add_variable("epsilon_xx",FIRST,LAGRANGE)); 
      _strain_variable_numbers.push_back(add_variable("epsilon_yy",FIRST,LAGRANGE)); 
      _strain_variable_numbers.push_back(add_variable("epsilon_zz",FIRST,LAGRANGE)); 
      _strain_variable_numbers.push_back(add_variable("epsilon_xy",FIRST,LAGRANGE)); 
      _strain_variable_numbers.push_back(add_variable("epsilon_xz",FIRST,LAGRANGE)); 
      _strain_variable_numbers.push_back(add_variable("epsilon_yz",FIRST,LAGRANGE)); 


      
      _strain_variables.resize(6); 
      _strain_variables[0]= "epsilon_xx"; 
      _strain_variables[1]= "epsilon_yy";
      _strain_variables[2] ="epsilon_zz"; 
      _strain_variables[3]= "epsilon_xy"; 
      _strain_variables[4]= "epsilon_xz";
      _strain_variables[5] ="epsilon_yz"; 
 
    }
  
}


void DiffCode::StressSystem::add_invariants(const std::string& invariant)
{
  _invariants.push_back(invariant); 
  if (!has_variable(invariant))
    {
      unsigned int var= add_variable(invariant,FIRST,LAGRANGE);
      libMesh::out<<var<<"\n";
    }
  


}


std::vector<std::string> DiffCode::StressSystem::get_invariants() 
{
  return _invariants; 
  
}


void DiffCode::StressSystem::compute_stresses_and_strains()
{
;
  // get the average to nodes system 
  
  AverageToNodes& a_to_n = *(_diff_code._average_to_nodes); 
  //  StressSolver& stress_solver = *(_diff_code._stress_solver); 

  
  a_to_n.compute_stresses_and_strains(
				      ); 
 
  stresses_computed = true;
}



Real DiffCode::StressSystem::get_strain_energy(DenseVector<Number>& stress_vector, 
					       DenseVector<Number>& strain_vector)
{
  Real strain_energy = 0 ; 

  for (unsigned int i=0; i< stress_vector.size(); i++) 
    strain_energy+=stress_vector(i)*strain_vector(i);
  
  return strain_energy/2.0; 

}


Real DiffCode::StressSystem::get_pressure(DenseVector<Number>& stress_vector) 
{
  Real pressure = 0; 
  Real n_stress = (Real) _diff_code._dim;
  for (unsigned int i=0; i< _diff_code._dim; i++) 
    pressure+=stress_vector(i); 
  
  pressure/=n_stress; 
  return -pressure; 
}


Real DiffCode::StressSystem::get_von_mises_stress(DenseVector<Number>& stress_vector) 
{

  Real von_mises_stress = 0; 
  if (_diff_code._dim == 3)
    {
      von_mises_stress = std::pow((stress_vector(0)-stress_vector(1)),2)+
	std::pow((stress_vector(1)-stress_vector(2)),2)+
	std::pow((stress_vector(2)-stress_vector(0)),2)+
	6*(std::pow(stress_vector(3),2) 
	   +std::pow(stress_vector(4),2)
	   +std::pow(stress_vector(5),2));
      von_mises_stress = std::sqrt(von_mises_stress/2);
    }
  else if (_diff_code._dim ==2) 
    {
      von_mises_stress =std::pow((stress_vector(0)-stress_vector(1)),2)+
	std::pow(stress_vector(0),2)+
	std::pow(stress_vector(1),2)+
	6*(std::pow(stress_vector(2),2));
      von_mises_stress =std::sqrt(von_mises_stress/2); 
    }

  return von_mises_stress;
}


void DiffCode::StressSystem::compute_invariants()
{

  if (!stresses_computed)
    compute_stresses_and_strains(); 
  

  std::vector<unsigned int> inv_variables;
  
  DenseVector<Number> strain_vector; 
  DenseVector<Number> stress_vector; 
  strain_vector.resize(3); 
  strain_vector.resize(3);

  if(_diff_code._dim==3)
    {
      strain_vector.resize(6); 
      stress_vector.resize(6);
    }


  for(unsigned int i=0; i<_invariants.size(); i++) 
    {
      libMesh::out<<"invariants "<<_invariants[i]<<"\n";
      
      const unsigned int var = variable_number(_invariants[i]);
      inv_variables.push_back(var);
    }

  unsigned int sys_num = number(); 
  MeshBase::const_node_iterator node_it = _diff_code._mesh->local_nodes_begin();
  const MeshBase::const_node_iterator node_it_end = _diff_code._mesh->local_nodes_end() ; 
  unsigned int n_calls = 0; 
  libMesh::out<<"aha"<<"\n";
  for(;node_it!=node_it_end;node_it++) 
    {
      Node* node =*node_it; 
      /// This check is done only once because the variables are all in the same places. 
      unsigned int n_comp = node->n_comp(sys_num,_stress_variable_numbers[0]) ;

      if (n_comp>0) 
	{
	  for (unsigned int i=0; i<_stress_variables.size(); i++) 
	    {
	      unsigned int dof_number = node->dof_number(sys_num,_stress_variable_numbers[i],0); 
	      // get the value from the solution vector 
	      stress_vector(i) = solution->el(dof_number); 
	    }
	  
	  for (unsigned int i=0; i<_strain_variables.size(); i++) 
	    {
	      unsigned int dof_number = node->dof_number(sys_num,_strain_variable_numbers[i],0); 
	      strain_vector(i) = solution->el(dof_number); 
	      
	    }
	  for (unsigned int i=0; i<inv_variables.size();i++)
	    {
	      Real ans =0; 
	      unsigned int dof_number = node->dof_number(sys_num,inv_variables[i],0);
	      if (_invariants[i]=="p")
		{
		  ans=get_pressure(stress_vector);
		
		}
	      else if (_invariants[i]=="ener") 
		{
		  
		  ans=get_strain_energy(stress_vector,strain_vector);
		  for (unsigned int i=0; i<_strain_variables.size();i++) 
		    {
		      libMesh::out<<_strain_variables[i]<<strain_vector(i)<<"\t";
		        libMesh::out<<_stress_variables[i]<<stress_vector(i)<<"\t";
		    }
		  libMesh::out<<"\n";

		  n_calls+=1;
		  libMesh::out<<"n "<<n_calls<<" ener "<<ans<<"\n";
		}
	      else if (_invariants[i]=="vm") 
		ans=get_von_mises_stress(stress_vector); 
	      else 
		ans=0; 
	      
	      
	      solution->set(dof_number,ans); 
	    }
	  
	}
    }
  libMesh::out<<"n_calls"<<n_calls<<"\n";


  // set stresses computed to false again
  stresses_computed=false;
}





