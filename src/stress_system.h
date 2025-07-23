/**
 * @file   stress_system.h
 * @author Subramanya Sadasiva <ssadasiv@me-98-16.dhcp.ecn.purdue.edu>
 * @date   Mon Mar 11 15:08:07 2013
 * 
 * @brief  This object holds the gradient field information at the nodes for all the 
 * stress driving forces in the system. For the non-stress systems, this information is stored in
 * a GradientSystem object. 
 */

#ifndef __stress_system_h__
#define __stress_system_h__


//libmesh includes 
#include "libmesh/explicit_system.h" 

#include "diff_code.h"
#include "average_to_nodes.h"
#include "stress_solver.h"
class DiffCode::StressSystem: public ExplicitSystem
{
 public:
  /** 
   * Constructor for the stress system. Used to store all the stress systems 
   * from the individual systems. 
   * @param eqSys  The parent equations systems object 
   * @param name  the name of this system 
   * @param number  the number of this system 
   * @return 
   */
  StressSystem(
		 EquationSystems& eqSys, 
		 const std::string& name, 
		 const unsigned int& number
	       );
  
  /** 
   * Function to add an invariant. 
   *  
   * @param invariant Invariant to be added  
   */
  void add_invariants(const std::string& invariant); 
 
  /** 
   * Compute the stresses and strains for the given stress system. 
   * 
   * @param stress_system 
   */
  void compute_stresses_and_strains();

  /** 
   * Compute the invariants for the stresses. 
   * param stress_system_number
   */
  void compute_invariants(); 
  
  /** 
   * 
   *  Return the vector containing the invariants. 
   *  I am going to compute the averaged von mises stress
   *  and the averaged strain energy in average to nodes itself. 
   * The trick is to try and make stress_system a friend class to the 
   * the other one. 
   * @return 
   */
  std::vector<std::string> get_invariants(); 


   /** 
   * Computes the first invariant of the stress given the 
   * stress vector 
   * @param stress_vector stress vector at the point. 
   * 
   * @return 
   */
  Real get_pressure(DenseVector<Number>& stress_vector);

  /** 
   * Computes the strain energy stored at the point 
   * 
   * @param stress_vector stress vector at the point 
   * @param strain_vector strain vector at the point
   * 
   * @return 
   */
  Real get_strain_energy(DenseVector<Number>& stress_vector, 
		       DenseVector<Number>& strain_vector); 


  /** 
   * Computes the von mises strain at the point . 
   * 
   * @param stress_vector stress vector at the point 
   * 
   * @return 
   */
  Real get_von_mises_stress(DenseVector<Number>& stress_vector); 


 private: 
  /// Reference to the parent diff_code object. 
  DiffCode& _diff_code; 
  
  /// system number of the average to nodes system. 
  //AverageToNodes&  _a_to_n;
  
  ///  Material - An object containing a stiffness matrix object. that can be computed 
  /// so that I can do the material properties easily. 
  std::vector<Material> _mats; 
  /// Vector containing all the invariants that need to be computed 
  std::vector<std::string> _invariants; 
  /// Vector containing the stress variables
  std::vector<std::string> _stress_variables; 
  /// Vector containing the strain variables
  std::vector<std::string> _strain_variables; 
  /// Vector containing the integers corresponding to the stress variables
  std::vector<unsigned int> _stress_variable_numbers; 
  /// Vector containing the integers corresponding to the strain variables 
  std::vector<unsigned int> _strain_variable_numbers; 


 
  
  /// Whether the stresses are computed 

  bool stresses_computed; 


  /// Stress system 
  /// StressSolver& _str_sys;
};

#endif 
