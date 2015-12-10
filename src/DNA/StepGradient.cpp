// StepGradient structure
// Nicolas Clauvelin


#include <IndexManager.h>
#include <StepGradient.h>


// flattening function for StepGradientVec
VectorN flatten_StepGradientVec(const StepGradientVec& v) {

    // stl container
    std::vector<Real> vec;
    vec.reserve(v.size()*emDNAConstants::StepParametersDim);

    // flattening
    std::for_each(v.begin(), v.end(), [&vec](const StepGradient& g){
        vec.insert(vec.end(),
                   g._rotation.packed_values()._values,
                   g._rotation.packed_values()._values+3);
        vec.insert(vec.end(),
                   g._translation.packed_values()._values,
                   g._translation.packed_values()._values+3);
    });

    // VectorN container
    return VectorN(vec);

};


// filtering function for a vector of StepGradient
VectorN filter_free_dofs_in_StepGradientVec(const StepGradientVec& v,
                                            const IndexManager& idx_mgr) {

    // flattening
    const VectorN vec = flatten_StepGradientVec(v);

    // free dofs indices
    const std::vector<Size> free_idx = idx_mgr.inline_free_dofs_indices();
    const Size n_idx = free_idx.size();

    // filtering
    // this assumes that the free dofs are sorted
    VectorN filtered_v(free_idx.size(), FLOAT_INIT);
    for (Size i=0; i<n_idx; ++i)
        filtered_v[i] = vec[free_idx[i]];

    return filtered_v;

};

