
// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// C++ includes
#include <string>
#include <cstdlib> // std::strtol
#include <sstream>

// Local includes
//#include "libmesh/abaqus_io.h"
#include "abaqus_io.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"
#include "libmesh/enum_elem_type.h"


using namespace libMesh;
// Anonymous namespace to hold mapping Data for Abaqus/libMesh element types
namespace
{
  /**
   * Data structure used for mapping Abaqus IDs to libMesh IDs, and
   * eventually (possibly) vice-versa.
   */
  struct ElementDefinition
  {
    // Maps (zero-based!) Abaqus local node numbers to libmesh local node numbers
    std::vector<unsigned> abaqus_zero_based_node_id_to_libmesh_node_id;

    // Maps (zero-based!) Abaqus side numbers to libmesh side numbers
    std::vector<unsigned> abaqus_zero_based_side_id_to_libmesh_side_id;
  };

  /**
   * Locally-available map containing all element data.
   */
  std::map<libMesh::ElemType, ElementDefinition> eletypes;

  /**
   * Helper function to fill up eletypes map
   */
  void add_eletype_entry(libMesh::ElemType libmesh_elem_type,
			 const unsigned* node_map,
			 unsigned node_map_size,
			 const unsigned* side_map,
			 unsigned side_map_size)
  {
    // If map entry does not exist, this will create it
    ElementDefinition& map_entry = eletypes[libmesh_elem_type];


    // Use the "swap trick" from Scott Meyer's "Effective STL" to swap
    // an unnamed temporary vector into the map_entry's vector.  Note:
    // the vector(iter, iter) constructor is used.
    std::vector<unsigned>(node_map,
			  node_map+node_map_size).swap(map_entry.abaqus_zero_based_node_id_to_libmesh_node_id);

    std::vector<unsigned>(side_map,
			  side_map+side_map_size).swap(map_entry.abaqus_zero_based_side_id_to_libmesh_side_id);
  }


  /**
   * Helper function to initialize the eletypes map.
   */
  void init_eletypes ()
  {
    // This should happen only once.  The first time this method is
    // called the eletypes data struture will be empty, and we will
    // fill it.  Any subsequent calls will find an initialized
    // eletypes map and will do nothing.
    if (eletypes.empty())
      {
	{
	  // TRI3
	  const unsigned int node_map[] = {0,1,2}; // identity
	  const unsigned int side_map[] = {0,1,2}; // identity
	  add_eletype_entry(TRI3, node_map, 3, side_map, 3);
	}

	{
	  // QUAD4
	  const unsigned int node_map[] = {0,1,2,3}; // identity
	  const unsigned int side_map[] = {0,1,2,3}; // identity
	  add_eletype_entry(QUAD4, node_map, 4, side_map, 4);
	}

	{
	  // TET4
	  const unsigned int node_map[] = {0,1,2,3}; // identity
	  const unsigned int side_map[] = {0,1,2,3}; // identity
	  add_eletype_entry(TET4, node_map, 4, side_map, 4);
	}

	{
	  // TET10
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9}; // identity
	  const unsigned int side_map[] = {0,1,2,3};             // identity
	  add_eletype_entry(TET10, node_map, 10, side_map, 4);
	}

	{
	  // HEX8
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7}; // identity
	  const unsigned int side_map[] = {0,5,1,2,3,4};     // inverse = 0,2,3,4,5,1
	  add_eletype_entry(HEX8, node_map, 8, side_map, 6);
	}

	{
	  // HEX20
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15}; // map is its own inverse
	  const unsigned int side_map[] = {0,5,1,2,3,4};                                       // inverse = 0,2,3,4,5,1
	  add_eletype_entry(HEX20, node_map, 20, side_map, 6);
	}

	{
	  // HEX27
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,26,20,25,21,22,23,24}; // inverse = ...,21,23,24,25,26,22,20
	  const unsigned int side_map[] = {0,5,1,2,3,4};                                                            // inverse = 0,2,3,4,5,1
	  add_eletype_entry(HEX27, node_map, 27, side_map, 6);
	}

	{
	  // PRISM6
	  const unsigned int node_map[] = {0,1,2,3,4,5}; // identity
	  const unsigned int side_map[] = {0,4,1,2,3};   // inverse = 0,2,3,4,1
	  add_eletype_entry(PRISM6, node_map, 6, side_map, 5);
	}

	{
	  // PRISM15
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11}; // map is its own inverse
	  const unsigned int side_map[] = {0,4,1,2,3};                          // inverse = 0,2,3,4,1
	  add_eletype_entry(PRISM15, node_map, 15, side_map, 5);
	}

	{
	  // PRISM18
	  const unsigned int node_map[] = {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11,15,16,17}; // map is its own inverse
	  const unsigned int side_map[] = {0,4,1,2,3};                                   // inverse = 0,2,3,4,1
	  add_eletype_entry(PRISM18, node_map, 18, side_map, 5);
	}



      } // if (eletypes.empty())
  } // init_eletypes()
} // anonymous namespace





namespace libMesh
{

  AbaqusIO::AbaqusIO (MeshBase& mesh_in) :
    MeshInput<MeshBase> (mesh_in),
    build_sidesets_from_nodesets(false),
    _already_seen_part(false)
  {
  }




  AbaqusIO::~AbaqusIO ()
  {
  }




  void AbaqusIO::read (const std::string& fname)
  {
    // Get a reference to the mesh we are reading
    MeshBase& the_mesh = MeshInput<MeshBase>::mesh();

    // Clear any existing mesh data
    the_mesh.clear();

    // Open stream for reading
    _in.open(fname.c_str());
    libmesh_assert(_in.good());

    // Read file line-by-line... this is based on a set of different
    // test input files.  I have not looked at the full input file
    // specs for Abaqus.
    std::string s;
    while (true)
      {
	// Try to read something.  This may set EOF!
	std::getline(_in, s);

	if (_in)
	  {
	    // Process s...
	    //
	    // There are many sections in Abaqus files, we read some
	    // but others are just ignored...  Some sections may occur
	    // more than once.  For example for a hybrid grid, you
	    // will have multiple *Element sections...

	    // Some Abaqus files use all upper-case for section names,
	    // so we will just convert s to uppercase
	    std::string upper(s);
	    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

	    // 0.) Look for the "*Part" section
	    if (upper.find("*PART") == 0)
	      {
		// libMesh::out << "Found parts section!" << std::endl;

		if (_already_seen_part)
		  {
		    libMesh::err << "We currently don't support reading Abaqus files with multiple PART sections" << std::endl;
		    libmesh_error();
		  }

		_already_seen_part = true;
	      }

	    // 1.) Look for the "*Nodes" section
	    if (upper.find("*NODE") == 0)
	      {
		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read a block of nodes
		this->read_nodes();
	      }



	    // 2.) Look for the "*Element" section
	    else if (upper.find("*ELEMENT,") == 0)
	      {
		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read a block of elements
		this->read_elements(upper);
	      }



	    // 3.) Look for a Nodeset section
	    else if (upper.find("*NSET") == 0)
	      {
		std::string nset_name = this->parse_label(s, "nset");
		
		bool generated = false; 
		// check if the nodeset is generated
		if (upper.find("GENERATE") !=std::string::npos) 
		  generated = true; 
		  
		
		// I haven't seen an unnamed elset yet, but let's detect it
		// just in case...
		if (nset_name == "")
		  {
		    libMesh::err << "Unnamed nset encountered!" << std::endl;
		    libmesh_error();
		  }
		
		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read the IDs, storing them in _nodeset_ids
		this->read_ids(nset_name, _nodeset_ids,generated);
	      } // *Nodeset



	    // 4.) Look for an Elset section
	    else if (upper.find("*ELSET") == 0)
	      {
		std::string elset_name = this->parse_label(s, "elset");
		// check if the element set is generated
		bool generated = false; 
		if (upper.find("GENERATE") !=std::string::npos) 
		  generated = true; 
	
		// internal subsets are not to be have there names set 
		// in the mesh. 
		bool set_subdomain_name = true;
		if (upper.find("INTERNAL")!=std::string::npos) 
		  set_subdomain_name = false; 
		
		_to_add_subdomain_id.insert(std::make_pair(elset_name,set_subdomain_name));
		
		// I haven't seen an unnamed elset yet, but let's detect it
		// just in case...
		if (elset_name == "")
		  {
		    libMesh::err << "Unnamed elset encountered!" << std::endl;
		    libmesh_error();
		  }

		// Debugging
		// libMesh::out << "Processing ELSET: " << elset_name << std::endl;

		// Process any lines of comments that may be present
		this->process_and_discard_comments();
		
		// Read the IDs, storing them in _elemset_ids
		this->read_ids(elset_name,_elemset_ids,generated);
		
	      } // *Elset



	    // 5.) Look for a Surface section.  Need to be a little
	    // careful, since there are also "surface interaction"
	    // sections we don't want to read here.
	    else if (upper.find("*SURFACE,") == 0)
	      {
		// libMesh::out << "Found SURFACE section: " << s << std::endl;
		
		// Get the name from the Name=Foo label.  This will be the map key.
		std::string sideset_name = this->parse_label(s, "name");
	     
		// Print name of section we just found
		// libMesh::out << "Found surface section named: " << sideset_name << std::endl;

		// Process any lines of comments that may be present
		this->process_and_discard_comments();

		// Read the sideset IDs
		this->read_sideset(sideset_name, _sideset_ids);

		// Debugging: print status of most recently read sideset
		// libMesh::out << "Read " << _sideset_ids[sideset_name].size() << " sides in " << sideset_name << std::endl;
	      }

	    continue;
	  } // if (_in)

	// If !file, check to see if EOF was set.  If so, break out
	// of while loop.
	if (_in.eof())
	  break;

	// If !in and !in.eof(), stream is in a bad state!
	libMesh::err << "Stream is bad!\n";
	libMesh::err << "Perhaps the file: " << fname << " does not exist?" << std::endl;
	libmesh_error();
      } // while


    //
    // We've read everything we can from the file at this point.  Now
    // do some more processing.
    //
    libMesh::out << "Mesh contains "
		 << the_mesh.n_elem()
		 << " elements, and "
		 << the_mesh.n_nodes()
		 << " nodes." << std::endl;

    // TODO: Remove these or write a function to do it?
//    {
//      container_t::iterator it=_nodeset_ids.begin();
//      for (; it != _nodeset_ids.end(); ++it)
//	{
//	  libMesh::out << "Node set '" << (*it).first << "' contains " << (*it).second.size() << " ID(s)." << std::endl;
//	}
//    }
//
//    {
//      container_t::iterator it=_elemset_ids.begin();
//      for (; it != _elemset_ids.end(); ++it)
//	{
//	  libMesh::out << "Elem set '" << (*it).first << "' contains " << (*it).second.size() << " ID(s)." << std::endl;
//	}
//    }


    // Set element IDs based on the element sets.
    this->assign_subdomain_ids();

    // Assign nodeset values to the BoundaryInfo object
    // this->assign_boundary_node_ids();

    // Assign sideset values in the BoundaryInfo object
    this->assign_sideset_ids();

    // Abaqus files only contain nodesets by default.  To be useful in
    // applying most types of BCs in libmesh, we will definitely need
    // sidesets.  So we can call the new BoundaryInfo function which
    // generates sidesets from nodesets.
    //if (build_sidesets_from_nodesets)
    // the_mesh.boundary_info->build_side_list_from_node_list();

  } // read()







  void AbaqusIO::read_nodes()
  {
    // Get a reference to the mesh we are reading
    MeshBase& the_mesh = MeshInput<MeshBase>::mesh();

    // Debugging: print node count
    // libMesh::out << "Before read_nodes(), mesh contains "
    // 		 << the_mesh.n_elem()
    // 		 << " elements, and "
    // 		 << the_mesh.n_nodes()
    // 		 << " nodes." << std::endl;

    // In the input file I have, Abaqus neither tells what
    // the mesh dimension is nor how many nodes it has...

    // The node line format is:
    // id, x, y, z
    // and you do have to parse out the commas.
    // The z-coordinate will only be present for 3D meshes

    // Temporary variables for parsing lines of text
    dof_id_type abaqus_node_id=0;
    Real x=0, y=0, z=0;
    char c;
    std::string dummy;

    // Defines the sequential node numbering used by libmesh
    dof_id_type libmesh_node_id = 0;

    // We will read nodes until the next line begins with *, since that will be the
    // next section.
    // TODO: Is Abaqus guaranteed to start the line with '*' or can there be leading white space?
    while (_in.peek() != '*' && _in.peek() != EOF)
      {
	// Re-Initialize variables to be read in from file
	abaqus_node_id=0;
	x = y = z = 0.;

	// Note: we assume *at least* 2D points here, should we worry about
	// trying to read 1D Abaqus meshes?
	_in >> abaqus_node_id >> c >> x >> c >> y;

	// Peek at the next character.  If it is a comma, then there is another
	// value to read!
	if (_in.peek() == ',')
	  _in >> c >> z;

	// Debugging: Print what we just read in.
	// libMesh::out << "node_id=" << node_id
	// 	     << ", x=" << x
	// 	     << ", y=" << y
	// 	     << ", z=" << z
	// 	     << std::endl;

	// Read (and discard) the rest of the line, including the newline.
	// This is required so that our 'peek()' at the beginning of this
	// loop doesn't read the newline character, for example.
	std::getline(_in, dummy);

	// Set up the abaqus -> libmesh node mapping.  This is usually just the
	// "off-by-one" map.
	_abaqus_to_libmesh_node_mapping[abaqus_node_id] = libmesh_node_id;

	// Add the point to the mesh using libmesh's numbering,
	// and post-increment the libmesh node counter.
	the_mesh.add_point(Point(x,y,z), libmesh_node_id++);
      } // while

    // Debugging: print node count.  Note: in serial mesh, this count may
    // be off by one, since Abaqus uses one-based numbering, and libmesh
    // just reports the length of its _nodes vector for the number of nodes.
    // libMesh::out << "After read_nodes(), mesh contains "
    //              << the_mesh.n_elem()
    //              << " elements, and "
    //              << the_mesh.n_nodes()
    //              << " nodes." << std::endl;

  } // read_nodes()





  void AbaqusIO::read_elements(std::string upper)
  {
    // Some *Element sections also specify an Elset name on the same line.
    // Look for one here.
    std::string elset_name = this->parse_label(upper, "elset");

    // Get a reference to the mesh we are reading
    MeshBase& the_mesh = MeshInput<MeshBase>::mesh();

    // initialize the eletypes map (eletypes is a file-global variable)
    init_eletypes();

    ElemType elem_type = INVALID_ELEM;
    unsigned n_nodes_per_elem = 0;

    // Within s, we should have "type=XXXX"
    if (upper.find("CPE4") != std::string::npos ||
	upper.find("CPS4") != std::string::npos)
      {
	elem_type = QUAD4;
	n_nodes_per_elem = 4;
	the_mesh.set_mesh_dimension(2);
      }
    else if (upper.find("CPS3") != std::string::npos)
      {
	elem_type = TRI3;
	n_nodes_per_elem = 3;
	the_mesh.set_mesh_dimension(2);
      }
    else if (upper.find("C3D8") != std::string::npos)
      {
	elem_type = HEX8;
	n_nodes_per_elem = 8;
	the_mesh.set_mesh_dimension(3);
      }
    else if (upper.find("C3D4") != std::string::npos)
      {
	elem_type = TET4;
	n_nodes_per_elem = 4;
	the_mesh.set_mesh_dimension(3);
      }
    else if (upper.find("C3D20") != std::string::npos)
      {
	elem_type = HEX20;
	n_nodes_per_elem = 20;
	the_mesh.set_mesh_dimension(3);
      }
    else if (upper.find("C3D6") != std::string::npos)
      {
	elem_type = PRISM6;
	n_nodes_per_elem = 6;
	the_mesh.set_mesh_dimension(3);
      }
    else if (upper.find("C3D15") != std::string::npos)
      {
	elem_type = PRISM15;
	n_nodes_per_elem = 15;
	the_mesh.set_mesh_dimension(3);
      }
    else if (upper.find("C3D10") != std::string::npos)
      {
	elem_type = TET10;
	n_nodes_per_elem = 10;
	the_mesh.set_mesh_dimension(3);
      }
    else
      {
	libMesh::err << "Unrecognized element type: " << upper << std::endl;
	libmesh_error();
      }

    // Insert the elem type we detected into the set of all elem types for this mesh
    _elem_types.insert(elem_type);

    // For reading in line endings
    std::string dummy;

    // Grab a reference to the element definition for this element type
    const ElementDefinition& eledef = eletypes[elem_type];

    // If the element definition was not found, the call above would have
    // created one with an uninitialized struct.  Check for that here...
    if (eledef.abaqus_zero_based_node_id_to_libmesh_node_id.size() == 0)
      {
	// libMesh::err << "No Abaqus->LibMesh mapping information for ElemType "
	// 	     << Utility::enum_to_string(elem_type)
	// 	     << "!"
	// 	     << std::endl;
	libmesh_error();
      }

    // We will read elements until the next line begins with *, since that will be the
    // next section.
    while (_in.peek() != '*' && _in.peek() != EOF)
      {
	// Read the element ID, it is the first number on each line.  It is
	// followed by a comma, so read that also.  We will need this ID later
	// when we try to assign subdomain IDs
	dof_id_type abaqus_elem_id = 0;
	char c;
	_in >> abaqus_elem_id >> c;

	// Debugging:
	// libMesh::out << "Reading data for element " << abaqus_elem_id << std::endl;

	// Add an element of the appropriate type to the Mesh.
	Elem* elem = the_mesh.add_elem(Elem::build(elem_type).release());

	// Associate the ID returned from libmesh with the abaqus element ID
	//_libmesh_to_abaqus_elem_mapping[elem->id()] = abaqus_elem_id;
	_abaqus_to_libmesh_elem_mapping[abaqus_elem_id] = elem->id();

	// The count of the total number of IDs read for the current element.
	unsigned id_count=0;

	// Continue reading line-by-line until we have read enough nodes for this element
	while (id_count < n_nodes_per_elem)
	  {
	    // Read entire line (up to carriage return) of comma-separated values
	    std::string csv_line;
	    std::getline(_in, csv_line);

	    // Create a stream object out of the current line
	    std::stringstream line_stream(csv_line);

	    // Process the comma-separated values
	    std::string cell;
	    while (std::getline(line_stream, cell, ','))
	      {
		// FIXME: factor out this strtol stuff into a utility function.
		char* endptr;
		long abaqus_global_node_id = std::strtol(cell.c_str(), &endptr, /*base=*/10);

		if (abaqus_global_node_id!=0 || cell.c_str() != endptr)
		  {
		    // Use the global node number mapping to determine the corresponding libmesh global node id
		    dof_id_type libmesh_global_node_id = _abaqus_to_libmesh_node_mapping[abaqus_global_node_id];

		    // Grab the node pointer from the mesh for this ID
		    Node* node = the_mesh.node_ptr(libmesh_global_node_id);

		    // Debugging:
		    // libMesh::out << "Assigning global node id: " << abaqus_global_node_id
		    //              << "(Abaqus), " << node->id() << "(LibMesh)" << std::endl;

		    // If node_ptr() returns NULL, it may mean we have not yet read the
		    // *Nodes section, though I assumed that always came before the *Elements section...
		    if (node == NULL)
		      {
			libMesh::err << "Error!  Mesh returned NULL Node pointer.\n";
			libMesh::err << "Either no node exists with ID " << libmesh_global_node_id
				     << " or perhaps this input file has *Elements defined before *Nodes?" << std::endl;
			libmesh_error();
		      }

		    // Note: id_count is the zero-based abaqus (elem local) node index.  We therefore map
		    // it to a libmesh elem local node index using the element definition map
		    unsigned libmesh_elem_local_node_id =
		      eledef.abaqus_zero_based_node_id_to_libmesh_node_id[id_count];

		    // Set this node pointer within the element.
		    elem->set_node(libmesh_elem_local_node_id) = node;

		    // Debugging:
		    // libMesh::out << "Setting elem " << elem->id()
		    //              << ", local node " << libmesh_elem_local_node_id
		    //              << " to global node " << node->id() << std::endl;

		    // Increment the count of IDs read for this element
		    id_count++;
		  } // end if strtol success
	      } // end while getline(',')
	  } // end while (id_count)

	// Ensure that we read *exactly* as many nodes as we were expecting to, no more.
	if (id_count != n_nodes_per_elem)
	  {
	    libMesh::err << "Error: Needed to read "
			 << n_nodes_per_elem
			 << " nodes, but read "
			 << id_count
			 << " instead!" << std::endl;
	    libmesh_error();
	  }

	// If we are recording Elset IDs, add this element to the correct set for later processing.
	// Make sure to add it with the Abaqus ID, not the libmesh one!
	if (elset_name != "")
	  {
	    // Debugging:
	    // libMesh::out << "Adding Elem " << abaqus_elem_id << " to Elmset " << elset_name << std::endl;
	    _elemset_ids[elset_name].push_back(abaqus_elem_id);
	  }
      } // end while (peek)
  } // read_elements()




  std::string AbaqusIO::parse_label(std::string line, std::string label_name)
  {
    // Do all string comparisons in upper-case
    std::string upper_line(line), upper_label_name(label_name);
    std::transform(upper_line.begin(), upper_line.end(), upper_line.begin(), ::toupper);
    std::transform(upper_label_name.begin(), upper_label_name.end(), upper_label_name.begin(), ::toupper);

    // Get index of start of "label="
    size_t label_index = upper_line.find(upper_label_name + "=");

    if (label_index != std::string::npos)
      {
	// Location of the first comma following "label="
	size_t comma_index = upper_line.find(",", label_index);

	// Construct iterators from which to build the sub-string.
	// Note the +1 is to skip past the "=" which follows the label name
	std::string::iterator
	  beg = line.begin() + label_name.size() + 1 + label_index,
	  end = (comma_index == std::string::npos) ? line.end() : line.begin()+comma_index;

	return std::string(beg, end);
      }

    // The label index was not found, return the empty string
    return std::string("");
  } // parse_label()


  void AbaqusIO::read_ids(std::string set_name, container_t& container, bool generated) 
  {
    
    if (generated) 
      {
	std::vector<dof_id_type>& id_storage  = container[set_name]; 
	while(_in.peek()!= '*' && _in.peek() !=EOF) 
	  {
	    // get the next line
	    std::string csv_line; 
	    std::getline(_in,csv_line); 
	   
	    // convert to sidestream
	    std::stringstream ss(csv_line); 
	    char c; 
	    // this is in the form init ',' final ',' step. 
	    dof_id_type init,final,step; 
	    ss>>init>>c>>final>>c>>step; 
	    for(dof_id_type i = init; i<final+1; i+=step) 
	      {
		id_storage.push_back(i); 
	      }
	  }
      }
    else
      {
	read_ids(set_name,container); 
      }
  }



  void AbaqusIO::read_ids(std::string set_name, container_t& container)
  {
    // Debugging
    // libMesh::out << "Reading ids for set: " << set_name << std::endl;

    // Grab a reference to a vector that will hold all the IDs
    std::vector<dof_id_type>& id_storage = container[set_name];

    // Read until the start of another section is detected, or EOF is encountered
    while (_in.peek() != '*' && _in.peek() != EOF)
      {
	// Read entire comma-separated line into a string
	std::string csv_line;
	std::getline(_in, csv_line);

	// On that line, use std::getline again to parse each
	// comma-separated entry.
	std::string cell;
	std::stringstream line_stream(csv_line);
	while (std::getline(line_stream, cell, ','))
	  {
	    // If no conversion can be performed by strtol, 0 is returned.
	    //
	    // If endptr is not NULL, strtol() stores the address of the
	    // first invalid character in *endptr.  If there were no
	    // digits at all, however, strtol() stores the original
	    // value of str in *endptr.
	    char* endptr;

	    // FIXME - this needs to be updated for 64-bit inputs
	    long id = std::strtol(cell.c_str(), &endptr, /*base=*/10);

	    // Note that lists of comma-separated values in abaqus also
	    // *end* with a comma, so the last call to getline on a given
	    // line will get an empty string, which we must detect.
	    if (id!=0 || cell.c_str() != endptr)
	      {
		// Debugging
		// libMesh::out << "Read id: " << id << std::endl;

		// 'cell' is now a string with an integer id in it
		id_storage.push_back( id );
	      }
	  }
      }

    // Status message
    // libMesh::out << "Read " << id_storage.size() << " ID(s) for the set " << set_name << std::endl;
  } // read_ids()




  void AbaqusIO::read_sideset(std::string sideset_name, sideset_container_t& container)
  {
    // Grab a reference to a vector that will hold all the IDs
    std::vector<std::pair<dof_id_type, unsigned> >& id_storage = container[sideset_name];

    // Variables for storing values read in from file
    dof_id_type elem_id=0;
    unsigned side_id=0;
    char c;
    
    // Read until the start of another section is detected, or EOF is encountered
    while (_in.peek() != '*' && _in.peek() != EOF)
      {
	// The strings are of the form: "391, S2" or "elset_name,S2" 
	std::string csv_line; 
	std::getline(_in,csv_line); 
	// get rid of the white space in the line .. there might be initial whitespaces
	csv_line.erase(std::remove_if(csv_line.begin(), 
				      csv_line.end(), 
				      ::isspace),
		       csv_line.end() ); 
	// get the first character in the line
	char first =csv_line.c_str()[0]; 
	
	// if the surface is specified using an element set, 
	// the first character is not a digit. 
	if (!::isdigit(first)) 
	  {
	    
	    std::stringstream ss(csv_line); 
	    std::vector<std::string> temp; 
	    // Tokenize into set_name and side definition
	    while (ss.good()) 
	      {
		std::string temp_string; 
		std::getline(ss,temp_string,','); 
		temp.push_back(temp_string); 
	      }
	    //get the name of the elset
	    std::string elset_name(temp[0]); 
	    // get the side definition from the file. 
	    std::stringstream ss2(temp[1]); 
	    
	    ss2>>c>>side_id; 
	    // get the element set corresponding to this sideset 
	    std::vector<dof_id_type>& id_vector = _elemset_ids[elset_name]; 
	    for(unsigned int i=0; i<id_vector.size(); i++) 
	      {
		
		id_storage.push_back(std::make_pair(id_vector[i],side_id)); 
	      }
	    
	  }
	else
	  {
	    // this is for the case where the surface is defined as an element number followed by the element "Elem_id,SX"
	    std::stringstream ss(csv_line); 
	    ss>>elem_id>>c; 
	    ss>>c>>side_id; 
	    id_storage.push_back(std::make_pair(elem_id,side_id)); 
	  } // while
      }
  }

    
    


  void AbaqusIO::assign_subdomain_ids()
  {
    // Get a reference to the mesh we are reading
    MeshBase& the_mesh = MeshInput<MeshBase>::mesh();

    // The number of elemsets we've found while reading
    std::size_t n_elemsets = _elemset_ids.size();


    unsigned int n_sub = 0; 


    // Fill in a temporary map with (ElemType, index) pairs based on the _elem_types set.  This
    // will allow us to easily look up this index in the loop below.
    std::map<ElemType, unsigned> elem_types_map;
    {
      unsigned ctr=0;
      for (std::set<ElemType>::iterator it=_elem_types.begin(); it!=_elem_types.end(); ++it)
	elem_types_map[*it] = ctr++;
    }

    // Loop over each Elemset and assign subdomain IDs to Mesh elements
    {
      // The elemset_id counter assigns a logical numbering to the _elemset_ids keys
      container_t::iterator it=_elemset_ids.begin();
      unsigned int elemset_id=0;
      for (;it!=_elemset_ids.end();++it)
	{
	  if (_to_add_subdomain_id[it->first]) 
	    {
	      elemset_id++;
	      /// add the things for the 
	      /// subdomain elemset by name map
	      

	      std::set<subdomain_id_type> subdomain_ids ; 
	      //for (std::set<ElemType>::iterator it2=_elem_types.begin();it2!=_elem_types.end();++it2)
	      //{
	      //  subdomain_ids.insert(elemset_id + (*it2)*n_elemsets);
	      //}
	     
	      
	      // Grab a reference to the vector of IDs
	      std::vector<dof_id_type>& id_vector = (*it).second;

	      // Loop over this vector
	      for (std::size_t i=0; i<id_vector.size(); ++i)
		{
		  // Map the id_vector[i]'th element ID (Abaqus numbering) to LibMesh numbering
		  dof_id_type libmesh_elem_id = _abaqus_to_libmesh_elem_mapping[ id_vector[i] ];
		  
		  // Get pointer to that element
		  Elem* elem = the_mesh.elem(libmesh_elem_id);
		  
		  if (elem == NULL)
		    {
		      libMesh::err << "Mesh returned NULL pointer for Elem " << libmesh_elem_id << std::endl;
		      libmesh_error();
		    }
		  // Assign this ID to the element in question
		  subdomain_id_type computed_id = elemset_id + (elem_types_map[elem->type()]*n_elemsets);
		  elem->subdomain_id() = computed_id;
		  subdomain_ids.insert(computed_id);
		  //	  _elemset_subdomain_map.insert(std::makepair(it->first,comuted_id));
		}
	      
	      _elemset_subdomain_map.insert(std::make_pair(it->first,subdomain_ids));
	    }
	}
    } 
  }// assign_subdomain_ids()



  void AbaqusIO::assign_boundary_node_ids()
  {
    // Get a reference to the mesh we are reading
    MeshBase& the_mesh = MeshInput<MeshBase>::mesh();

    // Iterate over the container of nodesets
    container_t::iterator it=_nodeset_ids.begin();
    for (unsigned current_id=0; it != _nodeset_ids.end(); ++it, ++current_id)
      {
	libMesh::out << "Assigning node boundary ID " << current_id << " to nodeset '"
		     << (*it).first
		     << "'." << std::endl;

	// Get a reference to the current vector of nodeset ID values
	std::vector<dof_id_type>& nodeset_ids = (*it).second;

	for (std::size_t i=0; i<nodeset_ids.size(); ++i)
	  {
	    // Map the Abaqus global node ID to the libmesh node ID
	    dof_id_type libmesh_global_node_id = _abaqus_to_libmesh_node_mapping[nodeset_ids[i]];

	    // Get node pointer from the mesh
	    Node* node = the_mesh.node_ptr(libmesh_global_node_id);

	    if (node == NULL)
	      {
		libMesh::err << "Error! Mesh returned NULL node pointer!" << std::endl;
		libmesh_error();
	      }

	    // Add this node with the current_id (which is determined by the
	    // alphabetical ordering of the map) to the BoundaryInfo object
	    the_mesh.boundary_info->add_node(node, current_id);
	  }
      }

  } // assign_boundary_node_ids()




  void AbaqusIO::assign_sideset_ids()
  {
    // Get a reference to the mesh we are reading
    MeshBase& the_mesh = MeshInput<MeshBase>::mesh();
   
    // initialize the eletypes map (eletypes is a file-global variable)
    init_eletypes();

    // Iterate over the container of sidesets
    sideset_container_t::iterator it=_sideset_ids.begin();
    for (unsigned current_id=0; it != _sideset_ids.end(); ++it, ++current_id)
      {
	libMesh::out << "Assigning sideset ID " << current_id << " to sideset '"
		     << (*it).first
		     << "'." << std::endl;

	// Get a reference to the current vector of nodeset ID values
	std::vector<std::pair<dof_id_type,unsigned> >& sideset_ids = (*it).second;

	for (std::size_t i=0; i<sideset_ids.size(); ++i)
	  {
	    // sideset_ids is a vector of pairs (elem id, side id).  Pull them out
	    // now to make the code below more readable.
	    dof_id_type  abaqus_elem_id = sideset_ids[i].first;
	    
	    unsigned abaqus_side_number = sideset_ids[i].second;

	    // Map the Abaqus element ID to LibMesh numbering
	    dof_id_type libmesh_elem_id = _abaqus_to_libmesh_elem_mapping[ abaqus_elem_id ];

	    // Get pointer to that element
	    Elem* elem = the_mesh.elem(libmesh_elem_id);


	    // Check that the pointer returned from the Mesh is non-NULL
	    if (elem == NULL)
	      {
		libMesh::err << "Mesh returned NULL pointer for Elem " << libmesh_elem_id << std::endl;
		libmesh_error();
	      }

	    // Grab a reference to the element definition for this element type
	    const ElementDefinition& eledef = eletypes[elem->type()];

	    // If the element definition was not found, the call above would have
	    // created one with an uninitialized struct.  Check for that here...
	    if (eledef.abaqus_zero_based_side_id_to_libmesh_side_id.size() == 0)
	      {
		libMesh::err << "No Abaqus->LibMesh mapping information for ElemType " << Utility::enum_to_string(elem->type()) << "!" << std::endl;
		libmesh_error();
	      }

	    // Add this node with the current_id (which is determined by the
	    // alphabetical ordering of the map).  Side numbers in Abaqus are 1-based,
	    // so we subtract 1 here before passing the abaqus side number to the
	    // mapping array
	    the_mesh.boundary_info->add_side(elem,
                                             eledef.abaqus_zero_based_side_id_to_libmesh_side_id[abaqus_side_number-1],
                                             current_id);

	    
	  }
	// set the name for the sideset. 

	the_mesh.boundary_info->sideset_name(current_id) = it->first; 
      }
  } // assign_sideset_ids()



  void AbaqusIO::process_and_discard_comments()
  {
    std::string dummy;
    while (true)
      {
	// We assume we are at the beginning of a line that may be
	// comments or may be data.  We need to only discard the line if
	// it begins with **, but we must avoid calling std::getline()
	// since there's no way to put that back.
	if (_in.peek() == '*')
	  {
	    // The first character was a star, so actually read it from the stream.
	    _in.get();

	    // Peek at the next character...
	    if (_in.peek() == '*')
	      {
		// OK, second character was star also, by definition this
		// line must be a comment!  Read the rest of the line and discard!
		std::getline(_in, dummy);

		// Debugging:
		// libMesh::out << "Read comment line: " << dummy << std::endl;
	      }
	    else
	      {
		// The second character was _not_ a star, so put back the first star
		// we pulled out so that the line can be parsed correctly by somebody
		// else!
		_in.unget();

		// Finally, break out of the while loop, we are done parsing comments
		break;
	      }
	  }
	else
	  {
	    // First character was not *, so this line must be data! Break out of the
	    // while loop!
	    break;
	  }
      }

  } // process_and_discard_comments()

 

  std::map<std::string,std::set<subdomain_id_type> > AbaqusIO::get_element_set_subdomain_ids()
  {
    return _elemset_subdomain_map; 
  }

} // namespace
