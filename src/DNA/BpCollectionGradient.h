// BpCollectionGradient functions
// Nicolas Clauvelin


// implementation of the matrices required for the computation of the gradient


#ifndef emDNA_BpCollectionGradient_h
#define emDNA_BpCollectionGradient_h


#include <StepGradient.h>
class BpStepParams;
class BpCollection;
class IndexManager;
class BpStepParamsGradient;


namespace BpCollectionGradient {


    // gradient computation functions
    VectorN free_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr);
    VectorN EEDR_collection_gradient(const BpCollection& bp_collection,
                                     const IndexManager& idx_mgr);
    VectorN EED_collection_gradient(const BpCollection& bp_collection,
                                    const IndexManager& idx_mgr);
    VectorN pulling_collection_gradient(const BpCollection& bp_collection,
                                        const IndexManager& idx_mgr,
                                        const Vector3& pulling_force);

    // single step parameters gradient function
    VectorN step_parameters_gradient(const BpStepParams& p,
                                     const BpStepParams& p0,
                                     const MatrixN& elastic_moduli);


    // following code is in an unnamed namespace
    namespace {


    // single step parameters gradient
    StepGradient single_step_parameters_gradient(const BpStepParams& p,
                                                 const BpStepParams& p0,
                                                 const MatrixN& elastic_moduli);

    // collection step parameters gradient
    // this gradient is computed over all free steps
    // (the step parameters gradient of a frozen step is zero by definition)
    StepGradientVec step_parameters_gradient(const BpCollection& bp_collection,
                                             const IndexManager& idx_mgr);

    // collection dofs gradient
    // this gradient is computed over all steps
    // this gradient is used as a basis for the computation of other gradients
    StepGradientVec dofs_gradient(const BpCollection& bp_collection,
                                  const IndexManager& idx_mgr);

    // imposed bp dofs gradient
    // this gradient is computed over all steps
    // the contribution of the step reduced by the boundary conditions is set to
    // zero
    StepGradientVec imposed_bp_dofs_gradient(const BpCollection& bp_collection,
                                             const IndexManager& idx_mgr);

    // imposed origin dofs gradient
    // imposed bp dofs gradient
    // this gradient is computed over all steps
    // the contribution of the step reduced by the boundary conditions is set to
    // zero (only for the translation component)
    StepGradientVec imposed_origin_dofs_gradient(const BpCollection&
                                                 bp_collection,
                                                 const IndexManager& idx_mgr);

    // frozen steps contributions
    StepGradientVec add_frozen_steps_contribs(const StepGradientVec& v,
                                              const BpCollection& bp_collection,
                                              const IndexManager& idx_mgr);


    }


}


#endif  // emDNA_BpCollectionGradient_h
