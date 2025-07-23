/**
 * @file   collated_phase_field.h
 * @author Subramanya <ssadasiv@me-98-2.dhcp.ecn.purdue.edu>
 * @date   Thu Oct 17 13:07:18 2013
 * 
 * @brief  This is a place to collate the results from all the 
 * phase field and diffusion systems that I have in my system. 
 * 
 * 
 */

#ifndef __collated_phase_field_h__
#define __collated_phase_field_h__
#include "libmesh/explicit_system.h"
#include "libmesh/id_types.h"
#include "diff_code.h" 

using namespace libMesh;
class DiffCode::CollatedPhaseField: public ExplicitSystem
{
 public: 
  CollatedPhaseField(
		     EquationSystems& eqSys, 
		     const std::string& name, 
		     const unsigned int number
		     ); 

  /** 
   * This is a function that takes the 
   * 
   */
  void collate_phase_field_values(); 

  /** 
   * This checks if the subdomain id given below is part of the subdomain ids over which the 
   * the phase field and the concentration field are defined. 
   * 
   * @param sub_id 
   * 
   * @return 
   */
  bool check_subdomain_inclusion(subdomain_id_type sub_id); 


  /** 
   * This gives me a material name if the subdomain name is matched against the 
   * 
   * 
   * @param sub_id 
   * 
   * @return 
   */
  std::string get_material_name(subdomain_id_type sub_id); 



  /** 
   * This is a function to refine the mesh based solely on the value of the field 
   * at the nodes. instead of the error vector. as the error vector based approach seems to be giving me
   * some fairly non-sensical results. 
   * let us see what happens because of this. 
   * 
   */
  void refine_mesh_based_on_values();

  void refine_mesh_based_on_error(); 


 private:
  DiffCode& _diff_code;   
  
  /// This is for subdomain ids so that it becomes easier to do things. 
  std::set<subdomain_id_type> _subdomain_ids; 
  /// I think I am going to add something here that makes it easy to get 
  /// the material properties. also I am going to assume that different 
  /// solder joints can have different material properties. 
  /// so I am going to add a map between the subdomain ids and the material 
  /// properties. so that I can check the things when I am solving the 
  /// the thing. 
  std::map<subdomain_id_type,std::string> _material_names; 

};




#endif 

