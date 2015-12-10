// BpCollectionGradient functions
// Nicolas Clauvelin


#include <IndexManager.h>
#include <BpCollection.h>
#include <BpGeometryFunctions.h>
#include <BpCollectionGeometry.h>
#include <BpCollectionGradient.h>


namespace BpCollectionGradient {


    // gradient computation function for free collection
    VectorN free_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr) {

        // dofs gradient
        const StepGradientVec grads = dofs_gradient(bp_collection, idx_mgr);

        // frozen step contributions
        const StepGradientVec frozen_grads =
        add_frozen_steps_contribs(grads, bp_collection, idx_mgr);

        // return the filtered gradient
        return filter_free_dofs_in_StepGradientVec(frozen_grads, idx_mgr);

    };


    // gradient computation function for collection with EEDR boundary
    // conditions
    VectorN EEDR_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr) {

        // dofs gradient
        const StepGradientVec grads = imposed_bp_dofs_gradient(bp_collection,
                                                               idx_mgr);

        // frozen step contributions
        const StepGradientVec frozen_grads =
        add_frozen_steps_contribs(grads, bp_collection, idx_mgr);

        // filtering
        VectorN grad_vec = filter_free_dofs_in_StepGradientVec(frozen_grads,
                                                               idx_mgr);

        // remove dofs related to boundary conditions
        // for EEDR boundary conditions this corresponds to the last six dofs
        return grad_vec.slice(0,
                              idx_mgr.n_of_free_dofs_variables()-
                              BoundaryConditionsDofsTrimming::EEDRCollection);

    };


    // gradient computation function for collection with EED boundary
    // conditions
    VectorN EED_collection_gradient(const BpCollection& bp_collection,
                                    const IndexManager& idx_mgr) {

        // dofs gradient
        const StepGradientVec grads =
        imposed_origin_dofs_gradient(bp_collection, idx_mgr);

        // frozen step contributions
        const StepGradientVec frozen_grads =
        add_frozen_steps_contribs(grads, bp_collection, idx_mgr);

        // filtering
        VectorN grad_vec = filter_free_dofs_in_StepGradientVec(frozen_grads,
                                                               idx_mgr);

        // remove dofs related to boundary conditions
        // for EED boundary conditions this corresponds to the last three dofs
        return grad_vec.slice(0,
                              idx_mgr.n_of_free_dofs_variables()-
                              BoundaryConditionsDofsTrimming::EEDCollection);

    };


    // gradient computation function for collection with pulling force on the
    // last bp
    VectorN pulling_collection_gradient(const BpCollection& bp_collection,
                                        const IndexManager& idx_mgr,
                                        const Vector3& pulling_force) {

        // dofs gradient
        StepGradientVec grads = dofs_gradient(bp_collection, idx_mgr);

        // add contributions from the pulling force
        for (Size i=0; i<bp_collection.n_of_bp_steps(); ++i) {
            grads[i]._translation[X] -=
            pulling_force[X]*emDNAConstants::PullingForceScaling;
            grads[i]._translation[Y] -=
            pulling_force[Y]*emDNAConstants::PullingForceScaling;
            grads[i]._translation[Z] -=
            pulling_force[Z]*emDNAConstants::PullingForceScaling;
        };

        // frozen steps contributions
        const StepGradientVec frozen_grads =
        add_frozen_steps_contribs(grads, bp_collection, idx_mgr);

        return filter_free_dofs_in_StepGradientVec(frozen_grads, idx_mgr);

    };


namespace {


    // single step parameters gradient
    StepGradient single_step_parameters_gradient(const BpStepParams& p,
                                                 const BpStepParams& p0,
                                                 const MatrixN& elastic_moduli)
    {

        // symmetrix elastic moduli matrix
        MatrixN symmetrize_fmat = elastic_moduli;
        symmetrize_fmat.transpose();
        symmetrize_fmat += elastic_moduli;

        // gradient
        VectorN g(Real(0.5)*(symmetrize_fmat*(p.inline_vector()-
                                              p0.inline_vector())));

        // gradient structure
        return StepGradient({
            { g[0]*RAD_2_DEG, g[1]*RAD_2_DEG, g[2]*RAD_2_DEG },
            { g[3], g[4], g[5] }
        });
        
    };


    // bp collection step parameters gradient
    // this vector is computed over all free steps
    StepGradientVec step_parameters_gradient(const BpCollection& bp_collection,
                                             const IndexManager& idx_mgr) {

        // free steps iterators
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();

        // container
        std::vector<StepGradient> grads;
        grads.reserve(idx_mgr.n_of_free_steps());

        // loop over free steps
        for (auto free_it = free_begin; free_it != free_end; ++free_it) {

            // step collection index
            const Size& coll_idx = free_it->collection_index();

            // step data
            const BpStepParams& p =
            bp_collection.bp_step_params(coll_idx);
            const BpStepParams& p0 =
            bp_collection.bp_step_intrinsic_parameters(coll_idx);
            const MatrixN& fmat =
            bp_collection.bp_step_force_constants(coll_idx);

            // single step gradient
            grads.push_back(std::
                            move(single_step_parameters_gradient(p, p0, fmat)));

        };

        return grads;

    };


    // collection dofs gradient
    // this gradient is computed over all steps
    StepGradientVec dofs_gradient(const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr) {

        // size
        const Size n_steps = bp_collection.n_of_bp_steps();

        // step parameters gradient
        // this gradient only contains free steps contributions
        const StepGradientVec sp_grad =
        std::move(step_parameters_gradient(bp_collection,
                                           idx_mgr));

        // free steps iterators
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();

        // dofs gradient container
        StepGradientVec dofs_grad(n_steps, StepGradient());

        // single step contributions
        // we loop over all free steps
        for (auto free_it = free_begin; free_it != free_end; ++free_it) {

            // collection index
            const Size& coll_idx = free_it->collection_index();

            // free step index
            const Size& free_idx = free_it->local_index();

            // matrices R and T
            const Matrix3 R_mat =
            BpCollectionGeometry::jacobian_matrix_R(coll_idx, bp_collection);
            const Matrix3 T_mat =
            BpCollectionGeometry::jacobian_matrix_T(coll_idx, bp_collection);

            // local terms
            // this corresponds to the single step dofs gradient
            StepGradient g = {
                sp_grad[free_idx]._rotation +
                    R_mat*sp_grad[free_idx]._translation,
                T_mat*sp_grad[free_idx]._translation
            };

            dofs_grad[coll_idx] = std::move(g);

        };

        // coupling contributions
        // we loop over all steps and check for free downward steps
        for (Size k=0; k<n_steps; ++k) {

            // we iterate over all free steps
            for (auto free_it = free_begin; free_it != free_end; ++free_it) {

                // collection index
                const Size& coll_idx = free_it->collection_index();

                // free step index
                const Size& free_idx = free_it->local_index();

                // check if this is a downward step
                // if yes, we compute the coupling contribution
                if (coll_idx > k) {

                    // matrix U
                    const Matrix3 U_mat =
                    BpCollectionGeometry::jacobian_matrix_U(k, coll_idx,
                                                            bp_collection);

                    // add the contribution to the gradient
                    dofs_grad[k]._rotation +=
                    U_mat*sp_grad[free_idx]._translation;

                };

            };

        };
		
		
		// add in electrostatic terms
		const float eps = 77.4; //permittivity of water at 300 K
		const float q = 1*(2*0.48*1.6021773e-19); //76% charge neutralized, 2 charges per step
		const float eps0 = 8.854188e-22;
		const float ion = 0.10; //10mL concentration, Debye length = 30.4A
		const float kap = 0.329*sqrt(ion); //Debye screening parameter
		//const float q = -0.24; //charge of DNA
		const float kT = 1.380658e-23*300.0;
		const float coeff = q*q/(4.0*M_PI*eps0*eps*kT);
		
		const std::vector<BpStepDofs>& bpdofs = bp_collection.bp_step_dofs();
		const std::vector<BasePair>& bps = bp_collection.base_pairs();
		BasePairConstIt begin = bps.begin();
		BasePairConstIt end = bps.end();
		
		/*//original code
		for (Size k=0; k<n_steps; ++k){
			// default constructor creates a zero vector
			for (auto i = begin; i<=begin+k; i++){				
				for (auto j=begin+k+1; j!= end; j++){
					Vector3 Dij = (*j).origin()-(*i).origin();
					Real r = Dij.norm();					
					dofs_grad[k]._translation += -coeff*exp(-kap*r)*(1/r)*(kap+1/r)*Dij*(1/r);
				}
			}
		}
		*/
		
		/*//week interaction		
		for (Size k=0; k<n_steps; k++){
			for (auto i = begin; i<begin+k; i=i+15){				
				for (auto j=end-((end-begin)%15); j> begin+k; j=j-15){
					Vector3 Dij = (*j).origin()-(*i).origin();
					Real r = Dij.norm();					
					if (r>30.395){//keep only the interaction beyond Debye length
						dofs_grad[k]._translation += -coeff*exp(-kap*r)*(1/r)*(kap+1/r)*Dij*(1/r);
					}
				}
			}
		}
		*/
		
		/*//omit interaction within Debye length
		for (Size k=0; k<n_steps; ++k){
			// default constructor creates a zero vector
			for (auto i = begin; i<=begin+k; i++){				
				for (auto j=begin+k+1; j!= end; j++){
					Vector3 Dij = (*j).origin()-(*i).origin();
					Real r = Dij.norm();
					if (r>30.395){//keep only the interaction beyond Debye length
						dofs_grad[k]._translation += -coeff*exp(-kap*r)*(1/r)*(kap+1/r)*Dij*(1/r);
					}
				}
			}
		}		
		*/
		
		
		//omit interaction within range
		for (Size k=0; k<n_steps; ++k){
			// default constructor creates a zero vector
			for (auto i = begin; i<=begin+k; i++){				
				for (auto j=begin+k+1; j!= end; j++){
					Vector3 Dij = (*j).origin()-(*i).origin();
					Real r = Dij.norm();
					if ((j-i)>=2){//keep only the interaction beyond range
						dofs_grad[k]._translation += -coeff*exp(-kap*r)*(1/r)*(kap+1/r)*Dij*(1/r);
					}
				}
			}
		}		
				
		//my code ends
		
        return dofs_grad;

    };


    // imposed bp dofs gradient
    StepGradientVec imposed_bp_dofs_gradient(const BpCollection& bp_collection,
                                             const IndexManager& idx_mgr) {

        // the reduced step is the last free step of the collection
        // we need the index manager to account for the case of the last step
        // being frozen
        const Size bc_step_idx = idx_mgr.bc_reduced_step().collection_index();

        // dofs gradient
        StepGradientVec dofs_grads = dofs_gradient(bp_collection, idx_mgr);

        // pre-computation of matrices K
        std::vector<Matrix3> K_matrices;
        K_matrices.reserve(bc_step_idx-1);
        for (Size i=0; i<bc_step_idx; ++i) {
            K_matrices.push_back(BpCollectionGeometry::
                                 jacobian_matrix_SigmaJ(bc_step_idx, i,
                                                        bp_collection));
            K_matrices.back().transpose();
        };

        // add boundary conditions contributions to the gradient
        for (Size i=0; i<bc_step_idx; ++i) {

            // boundary conditions contribution
            dofs_grads[i]._rotation -=
            K_matrices[i]*dofs_grads[bc_step_idx]._rotation;
            dofs_grads[i]._translation -=
            dofs_grads[bc_step_idx]._translation;

        };

        // the reduced step gradient is set to zero because of the boundary
        // conditions (imposed origin and orientation)
        dofs_grads[bc_step_idx] = StepGradient();

        return dofs_grads;

    };


    // imposed origin dofs gradient
    StepGradientVec imposed_origin_dofs_gradient(const BpCollection&
                                                 bp_collection,
                                                 const IndexManager& idx_mgr) {

        // the reduced step is the last free step of the collection
        // we need the index manager to account for the case of the last step
        // being frozen
        const Size bc_step_idx = idx_mgr.bc_reduced_step().collection_index();

        // dofs gradient
        StepGradientVec dofs_grads =
        dofs_gradient(bp_collection, idx_mgr);

        // add boundary conditions contributions to the gradient
        for (Size i=0; i<bc_step_idx; ++i) {

            // boundary conditions contribution
            dofs_grads[i]._translation -= dofs_grads[bc_step_idx]._translation;

        };

        // the reduced step gradient is set to zero because of the boundary
        // conditions (imposed origin and orientation)
        dofs_grads[bc_step_idx]._translation = Vector3();
        
        return dofs_grads;

    };


    // frozen steps contributions
    StepGradientVec add_frozen_steps_contribs(const StepGradientVec& grads,
                                              const BpCollection& bp_collection,
                                              const IndexManager& idx_mgr) {

        // nothing to do if no frozen steps
        if (idx_mgr.n_of_frozen_steps() == 0)
            return grads;

        // iterators
        const auto free_begin = idx_mgr.free_step_indexes_begin();
        const auto free_end = idx_mgr.free_step_indexes_end();
        const auto frozen_begin = idx_mgr.frozen_step_indexes_begin();
        const auto frozen_end = idx_mgr.frozen_step_indexes_end();

        // step gradients
        StepGradientVec step_grads = grads;

        // loop over all frozen steps
        for (auto frozen_it = frozen_begin; frozen_it != frozen_end;
             ++frozen_it) {

            // frozen step collection index
            const Size frozen_coll_idx = frozen_it->collection_index();

            // loop over all free steps located upward
            for (auto free_it = free_begin; free_it != free_end; ++free_it) {

                // check if we are dealing with an upward step
                if (free_it->collection_index() < frozen_coll_idx) {

                    // free step collection index
                    const Size free_coll_idx = free_it->collection_index();

                    // frozen jacobian (matrix W)
                    Matrix3 Wmat =
                    BpCollectionGeometry::
                    frozen_jacobian_matrix_W(free_coll_idx, frozen_coll_idx,
                                             bp_collection);
                    Wmat.transpose();

                    // add frozen terms
                    step_grads[free_coll_idx]._rotation +=
                    Wmat*grads[frozen_coll_idx]._translation;

                };

            };

            // set to zero the gradient of the frozen step
            step_grads[frozen_coll_idx] = StepGradient();
            
        };

        return step_grads;

    };


}


}

