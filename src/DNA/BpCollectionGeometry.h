// BpCollectionGeometry functions
// Nicolas Clauvelin


// functions for the geometry of bp collections
// see implementation notes for details


#ifndef emDNA_BpCollectionGeometry_h
#define emDNA_BpCollectionGeometry_h


#include <BpGeometryFunctions.h>
#include <BpCollection.h>
using namespace BpGeometryFunctions;


namespace BpCollectionGeometry {
    

    // ZYZ Euler angles from step parameters
    ZYZEulerAngles ZYZEulerAngles_from_step_parameters(const BpStepParams& p);

    // ZYZ Euler angles from step dofs
    ZYZEulerAngles ZYZEulerAngles_from_step_dofs(const BpStepDofs& dofs);

    // step frame
    Matrix3 step_frame_column_axes(Size step_index,
                                   const BpCollection& bp_collection);

    // jacobian matrix T
    Matrix3 jacobian_matrix_T(Size step_index,
                              const BpCollection& bp_collection);

    // jacobian matrix R
    Matrix3 jacobian_matrix_R(Size step_index,
                              const BpCollection& bp_collection);

    // infinitesimal rotation vectors S
    std::vector<Vector3> rotation_vectors_S(Size step_index,
                                            const BpCollection& bp_collection);

    // jacobian matrix U
    // index1 < index2
    Matrix3 jacobian_matrix_U(Size index1, Size index2,
                              const BpCollection& bp_collection);

    // jacobian matrix J
    // J^index1_index2
    // index1 > index2 (index1 has to be related to the bc reduced step)
    Matrix3 jacobian_matrix_J(Size index1, Size index2,
                              const BpCollection& bp_collection);

    // jacobian matrix SigmaJ
    // (Sigma^index1)^T.J^(index1)_(index2)
    // index1 > index2 (index1 has to be related to the bc reduced step)
    Matrix3 jacobian_matrix_SigmaJ(Size index1, Size index2,
                                   const BpCollection& bp_collection);

    // boundary conditions step matrix B for the different boundary conditions
    // bc_step_index is the collection index of the step reduced by the bc
    // this is a 6x6 matrix
    MatrixN bc_local_matrix_B_EEDR(Size bc_step_index, Size index,
                                   const BpCollection& bp_collection);
    MatrixN bc_local_matrix_B_EED(Size bc_step_index, Size index,
                                  const BpCollection& bp_collection);
    MatrixN bc_local_matrix_B_EER(Size bc_step_index, Size index,
                                  const BpCollection& bp_collection);

    // frozen step jacobian matrix W
    // index1 < index2
    Matrix3 frozen_jacobian_matrix_W(Size index1, Size index2,
                                     const BpCollection& bp_collection);
    MatrixN frozen_jacobian_matrix_C(Size index1, Size index2,
                                     const BpCollection& bp_collection);


}


#endif  // emDNA_BpCollectionGeometry_h
