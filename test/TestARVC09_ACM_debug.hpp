#include <cxxtest/TestSuite.h>

#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "HeartGeometryInformation.hpp"
#include "CardiacSimulationArchiver.hpp"
#include "HeartEventHandler.hpp"
#include "SingleTraceOutputModifier.hpp"
#include "RegularStimulus.hpp"
#include "ORdGksVarierARVC018.hpp"
#include "PostUpstrokeTimeAdaptivityControllerCL.hpp"
#include "Dijkstra_Cell_Factory_ACM_Endo.hpp"
#include "PseudoEcgCalculator.hpp"
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
#include "HeartConfig.hpp"
#include <Hdf5DataWriter.hpp>

#include "ConductivitiesModifierAnnulusIndexed3DExperimentalFibrosisSodium.hpp"
#include <algorithm> //for search of root node vec
#include <string>     // std::string, std::stoi

class TestSolveTorso : public CxxTest::TestSuite
{

public:
    void TestSolve() throw(Exception)
    {

	//Commandline arguments
   	std::vector <unsigned> nodes_present;
	unsigned dummy;
        bool root0 = CommandLineArguments::Instance()->OptionExists("--lv_basal_anterior_root");
        if (root0 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_basal_anterior_root");
            dummy = std::atoi(val);
			
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_basal_anterior_root not specified in commandline args" <<std::endl;
	  exit(1);
	}

        bool root1 = CommandLineArguments::Instance()->OptionExists("--lv_apical_posterior_root");
        if (root1 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_apical_posterior_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_apical_posterior_root not specified in commandline args" <<std::endl;
	  exit(1);
	}
        bool root2 = CommandLineArguments::Instance()->OptionExists("--lv_mid_posterior_root");
        if (root2 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_mid_posterior_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_mid_posterior_root not specified in commandline args" <<std::endl;
	  exit(1);
	}
        bool root3 = CommandLineArguments::Instance()->OptionExists("--lv_septal_root");
        if (root3 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_septal_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_septal_root not specified in commandline args" <<std::endl;
	  exit(1);
	}
        bool root4 = CommandLineArguments::Instance()->OptionExists("--rv_mid_anterolateral_root");
        if (root4 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--rv_mid_anterolateral_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "rv_mid_anterolateral_root not specified in commandline args, will place it anyways" <<std::endl;
	  nodes_present.push_back(1);

	}

        bool root5 = CommandLineArguments::Instance()->OptionExists("--rv_basal_posterolateral_root");
        if (root5 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--rv_basal_posterolateral_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "rv_basal_posterolateral_root not specified in commandline args" <<std::endl;
	  exit(1);
	}

        bool root6 = CommandLineArguments::Instance()->OptionExists("--rv_septal_root");
        if (root6 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--rv_septal_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "rv_septal_root not specified in commandline args" <<std::endl;
	  exit(1);
	}

       bool reg = CommandLineArguments::Instance()->OptionExists("--regions_file");
	std::string regions_file;
        if (reg == true) {
	
            regions_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--regions_file");

        }
	else
	{
	  std::cout << "regions_file not specified in commandline args" <<std::endl;
	  exit(1);
	}
	bool dummy_lv_root_nr = CommandLineArguments::Instance()->OptionExists("--nr_lv_roots");
	int nr_lv_roots;
        if (dummy_lv_root_nr == true) {
	
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--nr_lv_roots");
            nr_lv_roots = std::atoi(val);

        }
	else
	{
	  std::cout << "number of LV root nodes is not specified in commandline args" <<std::endl;
	  exit(1);
	}
		bool dummy_rv_root_nr = CommandLineArguments::Instance()->OptionExists("--nr_rv_roots");
	int nr_rv_roots;
        if (dummy_rv_root_nr == true) {
	
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--nr_rv_roots");
            nr_rv_roots = std::atoi(val);

        }
	else
	{
	  std::cout << "number of RV root nodes is not specified in commandline args" <<std::endl;
	  exit(1);
	}
	std::cout << "i get past commandline inputs" <<std::endl;
	std::string output_dir;
        output_dir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--save_dir");
        std::string apexbase_path;
        apexbase_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--apexbase");

        std::string mindists_path;
        mindists_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--mindists");




        std::string edge_nodes_path;
        edge_nodes_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--edge_nodes");
	    int fibrosis;
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--fibrosis");
            fibrosis = std::atoi(val);
           double scale_cond;
            char* val2 = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--scale_cond");
            scale_cond = std::atoi(val2);
           double time_step=0.008;
        //    char* val3 = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--pde_timestep");
          //  time_step = std::atoi(val3);



	/* Root node specification*/
		std::cout << "i get to root nodes" <<std::endl;
    	std::vector<unsigned int> rootnodes;    // a vector to hold apexbase data
	// open file    
	std::ifstream inputFile(edge_nodes_path+"roots/calibrated.txt");

	if (inputFile) {       
	int value;
	
	// read the elements in the file into a vector  
	while ( inputFile >> value ) {
            unsigned int value2 = value;
	rootnodes.push_back(value2);
	}
	}
  	std::vector<unsigned int> rootnodes_permuted;    // a vector to hold apexbase data
	// open file    
	std::ifstream inputFile1(edge_nodes_path+"roots/permuted.txt");
	inputFile.close();
	if (inputFile1) {       
	int value;
	
	// read the elements in the file into a vector  
	while ( inputFile1 >> value ) {
            unsigned int value2 = value;
	rootnodes_permuted.push_back(value2);
	}
	}
	inputFile1.close();
       // LV and RV activation root nodes
        std::vector<unsigned> lv_root_nodes;
	for (int i =0; i<nr_lv_roots;i++)
{
	if ( nodes_present[i]==1 )
	{
      	  lv_root_nodes.push_back(rootnodes[i]); // (5.577 9.896 6.501) LV anterior base
		  std::cout << "LV root " << lv_root_nodes[i] << std::endl;
	}
	  if( nodes_present[i] ==2)
	{
		 lv_root_nodes.push_back(rootnodes_permuted[i]); // (5.577 9.896 6.501) LV anterior base
	}  
	
}

        std::vector<unsigned> rv_root_nodes;
	for (int i =nr_lv_roots; i<nr_rv_roots;i++)
{
	if ( nodes_present[i]==1 )
	{
      	  rv_root_nodes.push_back(rootnodes[i]); // (5.577 9.896 6.501) LV anterior base
		  std::cout << "RV root " << rv_root_nodes[i] << std::endl;
	}
	  if( nodes_present[i] ==2)
	{
		 rv_root_nodes.push_back(rootnodes_permuted[i]); // (5.577 9.896 6.501) LV anterior base
	}  

}
	
    }
};
