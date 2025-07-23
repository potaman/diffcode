/*! New main file for the code. 
  
*/
//Include file 
#include "diff_code.h"
// Bring in everything from the libMesh namespace
using namespace libMesh;

// Print usage information if requested on command line
//void print_help(int argc, char** argv);

int main(int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);
  DiffCode* dc;
  DiffCode::Create(&dc,init.comm());
 
  dc->init();

  libMesh::out<<"test in the init"<<"\n";

  dc->run();

  DiffCode::Destroy(&dc);
  return 0;
}






