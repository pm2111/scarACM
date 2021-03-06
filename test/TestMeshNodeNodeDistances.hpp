

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"

#include "BidomainProblem.hpp"
#include "RegularStimulus.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCvodeCell.hpp"
#include "HeartEventHandler.hpp"
#include "HeartConfig.hpp"
#include "AbstractTetrahedralMesh.hpp"

#include "ORd2011epi_fkatpCvodeOpt.hpp"


#include "ConductivitiesModifierAnnulusIndexed.hpp"
#include <DistanceMapCalculator.hpp>


class TestConductivity : public CxxTest::TestSuite
{
public:  
    void TestConductivity2d() throw(Exception)
    {
		HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
		// MESH
		std::string filepath = "projects/pm2111/test/data/2Dmeshcm";
		//HeartConfig::Instance()->SetMeshFileName(filepath, cp::media_type::NoFibreOrientation);
                HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
		DistributedTetrahedralMesh<2,2> mesh;
		TrianglesMeshReader<2,2> mesh_reader(filepath);
                mesh.ConstructFromMeshReader(mesh_reader);
                HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
                HeartEventHandler::BeginEvent(HeartEventHandler::INITIALISE);

		const std::vector<unsigned> rNodePerm= mesh.rGetNodePermutation();
	        //create map which outputs the original node index, given the new permuted node index
	        std::cout <<" i get to the map creation bit" << std::endl;
		 std::map<unsigned,unsigned> MapNewOldIndex;
		 for( int i=0; i < rNodePerm.size() ; i++)
		{
			MapNewOldIndex[rNodePerm[i]] =i;
		}
			        std::cout <<" i get past the map creation bit" << std::endl;
		
		
        HeartEventHandler::EndEvent(HeartEventHandler::INITIALISE);
        // Problem immediately starts the EVERYTHING timer again
        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
		

	//Node to Node distances
	std::vector<double> node_distances; 
	DistanceMapCalculator<2,2> Calculator(mesh);
	std::vector<unsigned> target_node;
	target_node.push_back(8260); //-1 relative to matlab indexing
	Calculator.ComputeDistanceMap(target_node, node_distances);
	
    /* Save permutation just in case */
        OutputFileHandler output_file_handler("./annulus_indexed_fibrosis_new_stim/",false); // collective
        if ( PetscTools::AmMaster() )
        {
            out_stream perm_out_stream = output_file_handler.OpenOutputFile("/distances.txt", std::ios::out);
            for (unsigned i=0; i<node_distances.size(); ++i)
            {
                (*perm_out_stream) << node_distances[i] << "\n";
            }
        }
    



    }
};


