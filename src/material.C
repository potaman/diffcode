/**
 * @file   material.C
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Fri Apr 12 21:10:42 2013
 * 
 * @brief  Implementation of the material.h
 * 
 * I am implementing at present for the simplest kind of material. I am sure we can do this better.  
 */


#include "material.h" 
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"




Material::Material(const std::string& name, 
		     const unsigned int  number
		     )
{
  _name = name; 
  _number = number; 

}



RealTensor Material::get_diffusivity_tensor()
{
  // at present simple . later we will see if the temperature will need to be 
  // involved 

  return _diffusivity_tensor; 
}

DenseMatrix<Number>  Material::get_stiffness_matrix()
{
  return _stiffness_matrix;
}

RealTensor  Material::get_conductivity_tensor()
{
  return _conductivity_tensor; 
}


Real Material::get_surface_tension() 
{
  return _gamma; 
  
}

std::string Material::get_name() 
{
  return _name;
}

Real Material::get_conductivity() 
{
  return _conductivity;
  
}



Real Material::get_diffusivity()
{
  return _diffusivity;
}
  /** 
   * 
   * 
   * 
   * @return 
   */
  Real Material::get_atomic_volume()
  {
    return _atomic_volume; 
    
  }


  /** 
   * 
   * 
   * 
   * @return 
   */
  Real Material::get_surface_reaction_parameter()
  {
    return _surface_reaction_parameter; 
  }



  /** 
   * 
   * 
   * 
   * @return 
   */
  Real Material::get_surface_strain_energy_parameter()
  {
    return _surface_strain_energy_parameter; 
    
  }

  
  /** 
   * 
   * 
   * 
   * @return 
   */
  Real Material::get_electromigration_param()
  {
    return _electromigration_param; 
  }

  /** 
   * 
   * 
   * 
   * @return 
   */
  Real Material::get_surface_diffusivity()
  {
    return _surface_diffusivity; 

  }
  


  /** 
   * 
   *  
   * 
   * @return 
   */
  Real Material::get_surface_electromigration_param()
  {
    
    return _surface_electromigration_param; 
  }


//
//
void Material::get_stress(const std::vector<Real>& strain, 
			  std::vector<Real>& stress) 
{
  
  DenseVector<Number> stress_vec;
  stress_vec.resize(stress.size());  
  DenseVector<Number> strain_vec(strain); 
  
  _stiffness_matrix.vector_mult(stress_vec, strain_vec) ; 

  stress= stress_vec.get_values();  
}
//
void Material::get_thermal_strain(const Real& temperature, 
					    const Real& temp0,
					    std::vector<Real>& strain) 
{
  
  // Real alpha = parameters.get<Real>("alpha"); 
  Real strain_val = _thermal_expansivity*(temperature-temp0); 

  if (strain.size() == 3)
    {
      strain[0]= strain_val; 
      strain[1]= strain_val; 
      strain[2]= 0;  

    }
  else
    {
      strain[0]= strain_val; 
      strain[1]= strain_val; 
      strain[2]= strain_val;  
           
      strain[3]= 0.0; 
      strain[4]= 0.0; 
      strain[5]= 0.0;  

    }
}



void Material::get_concentration_strain(const std::vector<Real>& concentrations, 
						  const std::vector<Real>& conc_0,
						  std::vector<Real>& strain) 
{

  
  
  Real strain_val = _atomic_volume*(concentrations[0] - conc_0[0]) ;  
  
  
   if (strain.size() == 3)
    {
      strain[0]= strain_val; 
      strain[1]= strain_val; 
      strain[2]= 0;  

    }
  else
    {
      strain[0]= strain_val; 
      strain[1]= strain_val; 
      strain[2]= strain_val;  
           
      strain[3]= 0.0; 
      strain[4]= 0.0; 
      strain[5]= 0.0;  

    }
  



}


void Material::set_elasticity_parameters(Real youngs_modulus, 
					 Real poissons_ratio,
					 const int dim , 
					 bool plane_stress) 
{
  _youngs_modulus = youngs_modulus;
  _poissons_ratio = poissons_ratio;
  libMesh::out<<"poisson's ratio"<<_poissons_ratio;
  libMesh::out<<plane_stress<<std::endl;
  if (dim==2)
    {
      _stiffness_matrix.resize(3,3); 

      if (plane_stress) 
	{
	  Real c =  _youngs_modulus/(1-_poissons_ratio*_poissons_ratio);  
	  _stiffness_matrix(0,0) = 1.0; 
	  _stiffness_matrix(0,1) = _poissons_ratio; 
	  _stiffness_matrix(0,2) = 0.0;
 
	  _stiffness_matrix(1,0) = _poissons_ratio; 
	  _stiffness_matrix(1,1) = 1.0; 
	  _stiffness_matrix(1,2) = 0.0; 

	  _stiffness_matrix(2,0) = 0.0; 
	  _stiffness_matrix(2,1) = 0.0; 
	  _stiffness_matrix(2,2) = 1-_poissons_ratio/2.0; 

	  _stiffness_matrix*=c; 
	}
      else
	{
	  Real c =  _youngs_modulus/((1+_poissons_ratio)*(1-2*_poissons_ratio));
	  _stiffness_matrix(0,0) = 1-_poissons_ratio; 
	  _stiffness_matrix(0,1) = _poissons_ratio; 
	  _stiffness_matrix(0,2) = 0.0;
 
	  _stiffness_matrix(1,0) = _poissons_ratio; 
	  _stiffness_matrix(1,1) = 1-_poissons_ratio; 
	  _stiffness_matrix(1,2) = 0.0; 

	  _stiffness_matrix(2,0) = 0.0; 
	  _stiffness_matrix(2,1) = 0.0; 
	  _stiffness_matrix(2,2) = (1-2*_poissons_ratio)/2.0; 

	  _stiffness_matrix*=c; 

	}



    }
  else if (dim ==3) 
    {
      _stiffness_matrix.resize(6,6); 
      _stiffness_matrix*=0.0; 
      Real c =  _youngs_modulus/((1+_poissons_ratio)*(1-2*_poissons_ratio));

      _stiffness_matrix(0,0) = c*(1.0 - _poissons_ratio); 
      _stiffness_matrix(0,1) = c*(_poissons_ratio); 
      _stiffness_matrix(0,2) = c*(_poissons_ratio);
 
      _stiffness_matrix(1,0) = c*(_poissons_ratio); 
      _stiffness_matrix(1,1) = c*(1.0 - _poissons_ratio);
      _stiffness_matrix(1,2) = c*(_poissons_ratio); 

      _stiffness_matrix(2,0) = c*(_poissons_ratio); 
      _stiffness_matrix(2,1) = c*(_poissons_ratio); 
      _stiffness_matrix(2,2) = c*(1.0 - _poissons_ratio); 

      Real g = _youngs_modulus/2*(1+_poissons_ratio); 
      _stiffness_matrix(3,3) = g; 
      _stiffness_matrix(4,4) = g; 
      _stiffness_matrix(5,5) = g; 

    }
}


void Material::set_diffusivity(Real diffusivity)
{
  _diffusivity=diffusivity; 
}

void Material::set_diffusivity_tensor(RealTensor diffusivity_tensor)
{
  _diffusivity_tensor = diffusivity_tensor;
}

void Material::set_atomic_volume(Real atomic_volume) 
{
  _atomic_volume = atomic_volume;
}


void Material::set_thermal_conductivity(Real thermal_conductivity)
{
  _thermal_conductivity = thermal_conductivity; 
}


void Material::set_thermal_expansivity(Real thermal_expansivity) 
{
  _thermal_expansivity = thermal_expansivity ; 
}


void Material::set_conductivity(Real conductivity) 
{
  _conductivity = conductivity; 
} 


void  Material::set_surface_diffusivity(Real surface_diffusivity) 
{
  _surface_diffusivity = surface_diffusivity; 
}

void Material::set_electromigration_param(Real electromigration_param) 
{
  _electromigration_param = electromigration_param; 
}

void Material::set_thermomigration_param(Real thermomigration_param)
{
  _thermomigration_param = thermomigration_param; 
}

void Material::set_surface_thermomigration_param(Real surface_thermomigration_param)
{
  _surface_thermomigration_param = surface_thermomigration_param;
}

void Material::set_surface_reaction_param(Real surface_reaction_param) 
{
  _surface_reaction_parameter = surface_reaction_param; 
}

void Material::set_surface_strain_energy_param(Real surface_strain_energy_parameter)
{
  _surface_strain_energy_parameter = surface_strain_energy_parameter; 
}

void Material::set_surface_electromigration_param(Real surface_electromigration_param)
{
  _surface_electromigration_param = surface_electromigration_param;

}

void Material::set_surface_tension(Real surface_tension)
{
  _gamma = surface_tension;
  
}

void Material::set_heat_factor(Real heat_factor)
{
  _heat_factor = heat_factor; 
}

void Material::set_thermal_mass(Real thermal_mass)
{
  _thermal_mass = thermal_mass;   
}

Real Material::get_thermal_mass()
{
  return _thermal_mass; 
}

Real Material::get_thermal_conductivity()
{
  return _thermal_conductivity; 
}

Real Material::get_thermal_expansivity()
{
  return _thermal_expansivity; 
}

Real Material::get_heat_factor()
{
  return _heat_factor; 
}


Real Material::get_thermomigration_param()
{
  return _thermomigration_param; 
}

Real Material::get_surface_thermomigration_param()
{

  return _surface_thermomigration_param; 
}
