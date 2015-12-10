// BpCollectionElectrostaticEnergy class
// Juan Wei
#include <BpCollection.h>
#include <BpCollectionElectrostaticEnergy.h>


  const float eps = 77.4; //permittivity of water at 300 K
  const float q = 1*(2*0.48*1.6021773e-19); //76% charge neutralized, 2 charges per step
  const float eps0 = 8.854188e-22;
  const float ion = 0.10; //10mL concentration, Debye length = 30.4A
  const float kap = 0.329*sqrt(ion); //Debye screening parameter
  //const float q = -0.24; //charge of DNA
  const float kT = 1.380658e-23*300.0;
  const float coeff = q*q/(4.0*M_PI*eps0*eps*kT);
  
Real
BpCollectionElectrostaticEnergy::electrostatic_energy(const BpCollection& bp_collection) {
    
    Real E(FLOAT_INIT);
	
	const std::vector<BasePair>& bps = bp_collection.base_pairs();
	BasePairConstIt begin = bps.begin();
	BasePairConstIt end = bps.end();

	/*//original code
	for (auto i = begin; i< end; i++) { //for (auto i = begin; i< end-2; i++) (auto i = begin; i< end; i++)

		//BasePair& bpi = *i;
		
		for (auto j = i+1; j < end; j++) { //j=end causes energy infinity for (auto j = i+2; j < end; j++) (auto j = i+1; j < end; j++)
			//BasePair& bpj = *j;
			Vector3 Dij = (*j).origin()-(*i).origin();
			Real r = Dij.norm();	
			E += exp(-kap*r)/r; 
		};
		
		
    };    
	*/
	
	/*//weak interaction
	for (auto i = begin; i< end-15; i=i+15) { //for (auto i = begin; i< end-2; i++)		
		for (auto j = i+15; j < end; j=j+15) { //j=end causes energy infinity for (auto j = i+2; j < end; j++)
			Vector3 Dij = (*j).origin()-(*i).origin();
			Real r = Dij.norm();	
			if (r>30.395){//keep only the interaction beyond Debye length
				E += exp(-kap*r)/r; 
			}
		};
		
		
    };
	*/
	
	/*//omit interaction within Debye length
	for (auto i = begin; i< end; i++) { //for (auto i = begin; i< end-2; i++)
		for (auto j = i+1; j < end; j++) { //j=end causes energy infinity for (auto j = i+2; j < end; j++)
			Vector3 Dij = (*j).origin()-(*i).origin();
			Real r = Dij.norm();
			if (r>30.395){//keep only the interaction beyond Debye length
				E += exp(-kap*r)/r; 
			}
		}
		
		
    } 	
	*/
	
	//omit interaction within range
	for (auto i = begin; i< end; i++) { //for (auto i = begin; i< end-2; i++)
		for (auto j = i+2; j < end; j++) { //j=end causes energy infinity for (auto j = i+2; j < end; j++)
			Vector3 Dij = (*j).origin()-(*i).origin();
			Real r = Dij.norm();			
			E += exp(-kap*r)/r; 			
		}
		
		
    } 	
		
	E *= coeff; 

    return E;

};
