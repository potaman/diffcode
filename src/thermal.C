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


#include "thermal.h"
#include "collated_phase_field.h"

DiffCode::Thermal::Thermal(EquationSystems& eqSys,
				 const std::string& name, 
				 const unsigned int number):
  TransientLinearImplicitSystem(eqSys, name, number), 
  _diff_code(dynamic_cast<DiffCode&>(eqSys))
{
  const unsigned int dim = _diff_code._dim; 
  add_variable("t",FIRST,LAGRANGE); 
  _constrained=false; 
  attach_assemble_object(*this);
  attach_init_object(*this);

}

void DiffCode::Thermal::add_dirichlet_boundary_condition(boundary_id_type& i, 
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
  variables[0] = variable_number("t"); 

  // Create the boundary condition. 
  DirichletBoundary* dbc = new DirichletBoundary(bc_sides, 
						 variables, 
						 bc_func) ; 

  //Insert the boundary condition ito the set of boundary condiiton. 
  // this does not seem to be needed. need to figure out what to do . 
  _thermal_dirichlet_bcs.insert(std::make_pair(i,*dbc));
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
void DiffCode::Thermal::add_neumann_boundary_condition(boundary_id_type& i, 
							      RealGradient vector_value,
							      Real value ) 
  {
    // Insert the boundary condition sides into the neumann boundary condition set. 
    _dirichlet_boundary_sides.insert(i); 
    // Create the funcition for the boundary condition,
    FunctionBase<Real>* bc_func = new ConstFunction<Real>(value);
    // Create the function for the boundary condition
    FunctionBase<RealGradient>* bc_vec_func = new ConstFunction<RealGradient>(vector_value);   
    // Create a Neumann Boundary object
    NeumannBoundary* neumann_bc = new NeumannBoundary(i, 
						      *bc_vec_func,
						      *bc_func) ; 
    
    // Insert this into the map of boundary conditions. 
    _thermal_neumann_bcs.insert(std::make_pair(i,*neumann_bc));
    libMesh::out<<"Added boundary condition"<<"\n";
  }

//=================================================================================
// End of the neumann boundary condition
//=================================================================================

void DiffCode::Thermal::add_convection_boundary_condition(boundary_id_type& i, 
							    Real co_eff,
							    Real ambient_temp ) 
{
    // Insert the boundary condition sides into the neumann boundary condition set. 
  _convection_boundary_sides.insert(i); 
  ConvectionBC* cbc = new ConvectionBC();
  cbc->bc_side=i;
  cbc->ambient_temp=ambient_temp;
  cbc->co_eff = co_eff;
  // Insert this into the map of boundary conditions. 
  _thermal_convection_bcs.insert(std::make_pair(i,*cbc));
  libMesh::out<<"Added  convection boundary condition"<<"\n";
}
					  


//=================================================================================
// Initializer for the linear solver. 
//=================================================================================

void DiffCode::Thermal::reinit_lin_solver()
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
  ierr =KSPSetOptionsPrefix(solver_ksp,"temp_"); 
  
  // Set options from the options file 
  ierr =KSPSetFromOptions(solver_ksp); 
  
  // Put the solver back into the nonlinear solver 
  //linear_solver=petsc_solver; 

  release_linear_solver(petsc_solver); 
  
}


void DiffCode::Thermal::initialize()
{
  Parameters parameters; 
  parameters.set<Real>("t0") =_diff_code.get_equilibrium_temperature();
  libMesh::out<<"t0"<<_diff_code.get_equilibrium_temperature();
  project_solution(
		   DiffCode::Thermal::initial_temp,
		   NULL,
		   parameters); 
  *(old_local_solution)= *(current_local_solution); 
}

Number DiffCode::Thermal::initial_temp(const Point& p,
				    const Parameters& parameters,
				    const std::string&,
				    const std::string&)
{
  return parameters.get<Real>("t0"); 
}

void DiffCode::Thermal::assemble()
{

  Number dt = _diff_code._dt;

 const unsigned int dim = _diff_code._dim; 
  
  Real small_number = _diff_code._small_number; 
  CollatedPhaseField& ch = *(_diff_code._collated_phase_field);
    
  const NumericVector<Number>& ch_vector = *(ch.current_local_solution); 

  const ExplicitSystem& es = _diff_code.get_system<ExplicitSystem>(std::string("DiffCode::elec"));
  const NumericVector<Number>& e_vector = *(es.current_local_solution);


  const unsigned int t_var = variable_number("t");
  const unsigned int e_var = es.variable_number("e"); 
  const unsigned int phi_var = ch.variable_number("phi"); 
  
  const DofMap& dof_map_t = get_dof_map(); 
  std::vector<unsigned int> dof_indices_t; 

  const DofMap& dof_map_e = es.get_dof_map(); 
  std::vector<unsigned int> dof_indices_e; 
    
  const DofMap& dof_map_ch = ch.get_dof_map(); 
  std::vector<unsigned int> dof_indices_ch_phi; 

  FEType fe_e = es.variable_type(e_var); 
  FEType fe_t = variable_type(t_var);
  FEType fe_ch = ch.variable_type(phi_var); 

  AutoPtr<FEBase> fe_base_t(FEBase::build(dim,fe_t));
  AutoPtr<FEBase> fe_base_e(FEBase::build(dim,fe_e)); 
  AutoPtr<FEBase> fe_base_ch(FEBase::build(dim,fe_ch)); 

  QGauss q_bulk(dim, fe_e.default_quadrature_order()); 
  
  fe_base_t->attach_quadrature_rule(&q_bulk);
  fe_base_e->attach_quadrature_rule(&q_bulk); 
  fe_base_ch->attach_quadrature_rule(&q_bulk); 

  const std::vector<Real>& JxW = fe_base_e->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_ch = fe_base_ch->get_phi();  
  const std::vector<std::vector<Real> >& phi_e = fe_base_e->get_phi();
  const std::vector<std::vector<Real> >& phi_t = fe_base_t->get_phi(); 

  const std::vector<std::vector<RealGradient> >& dphi_e = fe_base_e->get_dphi(); 
  const std::vector<std::vector<RealGradient> >& dphi_t = fe_base_t->get_dphi(); 
  
  
  AutoPtr<FEBase> fe_base_surf_t(FEBase::build(dim,fe_t)); 
  QGauss q_face(dim-1,fe_t.default_quadrature_order()); 
 
  fe_base_surf_t->attach_quadrature_rule(&q_face); 
  
  const std::vector<Real>& JxW_surf = fe_base_surf_t->get_JxW(); 
  const std::vector<std::vector<Real> >& phi_surf = fe_base_surf_t->get_phi(); 
  const std::vector<Point>& xyz_surf = fe_base_surf_t->get_xyz(); 
  const std::vector<Point>& norm_surf = fe_base_surf_t->get_normals(); 
  
  
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
      
      dof_map_t.dof_indices(elem,dof_indices_t) ;
      dof_map_e.dof_indices(elem,dof_indices_e) ; 
      
      const unsigned int n_dofs = dof_indices_t.size(); 
      const unsigned int n_dofs_e = dof_indices_e.size(); 
      
      
      unsigned int n_dofs_ch; 
            
      if (ch.check_subdomain_inclusion(sub_id))
	{
	  dof_map_ch.dof_indices(elem, dof_indices_ch_phi, phi_var); 
	  n_dofs_ch = dof_indices_ch_phi.size(); 
	}
     
      Fe.resize(n_dofs); 
      Ke.resize(n_dofs,n_dofs); 

      fe_base_t->reinit(elem); 
      fe_base_e->reinit(elem);  
      fe_base_ch->reinit(elem); 
      
      // compute the value of the phase field variables
      std::string mat_name = ch.get_material_name(sub_id); 

      mat_it = _diff_code._materials.find(mat_name); 
      Real t_cond = mat_it->second.get_thermal_conductivity();
      Real t_mass = mat_it->second.get_thermal_mass(); 
      Real cond   = mat_it->second.get_conductivity(); 
      Real hf   = mat_it->second.get_heat_factor(); 


      for (unsigned int qp=0; qp<q_bulk.n_points(); qp++)
 
	{
	  // If you are in the active subdomain. compute the phi_qp and then compute the conductivity
	  Number t_old =0.0;
	  // this is the term for the electric potential gradient
	  Gradient grad_e; 
	  for (unsigned int i=0; i<n_dofs;i++)
	    {
	      t_old +=phi_t[i][qp]*old_local_solution->el(dof_indices_t[i]);
	    }
	  for (unsigned int i=0; i<n_dofs_e;i++)
	    {
	      grad_e.add_scaled(dphi_e[i][qp],e_vector(dof_indices_e[i]));
	    }
	  
	  if(ch.check_subdomain_inclusion(sub_id))
	    {
	      Real phi_qp_r = 0.0; 
	      for (unsigned int i = 0;i< n_dofs_ch;i++) 
		phi_qp_r += phi_ch[i][qp]*ch_vector(dof_indices_ch_phi[i]); 
	      Real heaviside = std::max((1+phi_qp_r)/2.0,0.0) +small_number ; 
	      t_cond*=heaviside;
	      cond*= std::max((1+phi_qp_r)/2.0,0.0); 
	      hf*= std::max((1+phi_qp_r)/2.0,0.0);
	      t_mass*=std::max((1+phi_qp_r)/2.0,0.0);
	    }
	  //	  libMesh::out<<hf<<std::endl;
	  //Compute the matrix element
	  for (unsigned int i=0; i< n_dofs;i++) 
	    {
	      Fe(i) +=JxW[qp]*(
			       t_mass*phi_t[i][qp]*t_old
			       +dt*phi_t[i][qp]*hf*cond*grad_e.size_sq()
			       );
	      
	      for(unsigned int j=0; j< n_dofs; j++) 
		Ke(i,j)+= JxW[qp]*( 
				   t_mass*phi_t[i][qp]*phi_t[j][qp]
				    +t_cond*dt*dphi_t[i][qp]*dphi_t[j][qp] 
				    );
	      
	    }

	}
      for (unsigned int side=0; side<elem->n_sides(); side++) 
	{
	  if (elem->neighbor(side) ==NULL) 
	    {
	      std::map<boundary_id_type,NeumannBoundary>::iterator it_nbc;
	      for (it_nbc=_thermal_neumann_bcs.begin();
		   it_nbc!=_thermal_neumann_bcs.end();
		   ++it_nbc)
		{
		  if (_diff_code._mesh->boundary_info->has_boundary_id(elem,side,it_nbc->first))
		    {
		      
		      NeumannBoundary* nbc = &(it_nbc->second);
		      if(nbc->is_active())
			{
			  fe_base_surf_t->reinit(elem,side); 
			  for(unsigned int qp =0; qp<q_face.n_points(); qp++) 
			    {
			      Real val = nbc->get_boundary_condition_value(xyz_surf[qp],
									   norm_surf[qp]);
			      for (unsigned int i=0;i<n_dofs;i++) 
				Fe(i)+=JxW_surf[qp]*val*phi_surf[i][qp];
			    }
			}
		    }
		}
	     
	      std::map<boundary_id_type,ConvectionBC>::iterator it_cbc;
	      for (it_cbc=_thermal_convection_bcs.begin();
	      	   it_cbc!=_thermal_convection_bcs.end();
	      	   ++it_cbc)
	      	{
	      	  if  (_diff_code._mesh->boundary_info->has_boundary_id(elem,side,it_cbc->first))
       
	      	    {
	          
	      	      //    libMesh::out<<"getting here at all"<<it_cbc->first;
	      	      fe_base_surf_t->reinit(elem,side); 
	      	      Number h = it_cbc->second.co_eff; 
	      	      Number at = it_cbc->second.ambient_temp; 
	      	      for (unsigned int qp = 0; qp<q_face.n_points(); qp++) 
	      		{
	     	  
	      		  for (unsigned int i= 0;i<n_dofs;i++)
	      		    {
	      		      Fe(i)+=JxW_surf[qp]*dt*phi_surf[i][qp]*h*at;
	      		      for (unsigned int j =0 ; j<n_dofs;j++) 
	      			Ke(i,j)+=JxW_surf[qp]*dt*phi_surf[i][qp]*phi_surf[j][qp]*h;
	     	      
	      		    }
	      		}
	      	    }
	      	}
	    }
	}
      dof_map_t.heterogenously_constrain_element_matrix_and_vector(Ke,Fe,dof_indices_t); 
      matrix->add_matrix(Ke,dof_indices_t); 
      rhs->add_vector(Fe,dof_indices_t); 
    }
  matrix->close(); 
  rhs->close(); 
  libMesh::out<<"assembly complete"<<std::endl;
}



void DiffCode::Thermal::solve_and_constrain()
{
  solve(); 
  get_dof_map().enforce_constraints_exactly(*this,solution.get());
  update();
  
}
