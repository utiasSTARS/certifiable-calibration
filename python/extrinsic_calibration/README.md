# Globally Optimal Hand-Eye Calibration
Certifiably globally optimal extrinsic calibration of a monocular camera to a sensor that provides egomotion estimates.

<img src="https://github.com/utiasSTARS/global-hand-eye/blob/master/python/extrinsic_calibration/egomotion.png" width="500">

## Installation and Dependencies

All code was developed in Python 3.5.

The following packages were used to compile and run the code:
- matplotlib~=3.0.3
- scipy~=1.2.1
- numpy~=1.16.3
- [liegroups~=1.0.0](https://github.com/utiasSTARS/liegroups.git)
- future~=0.18.2
- [cvxpy~=1.0.29](https://www.cvxpy.org/)

The [CVXOPT](http://cvxopt.org/) solver was used in this package.

## Minimal Usage
The script `example_hand_eye.py` demonstrates a sample usage of the solver on simulated noisy data.

## Repeating Experiments
The script `synthetic_experiments.py` runs the dual solver experiments shown in the paper.

## Citation
If you use any of this code in your work, please cite the [relevant publication](https://arxiv.org/abs/2005.08298):

```bibtex
@misc{wise2020certifiably,
    title={Certifiably Optimal Monocular Hand-Eye Calibration},
    author={Emmett Wise and Matthew Giamou and Soroush Khoubyarian and Abhinav Grover and Jonathan Kelly},
    year={2020},
    eprint={2005.08298},
    archivePrefix={arXiv},
    primaryClass={cs.RO}
}

```

## File Naming Convention
The result file names follow the convention:

`datasetfilename_constraints.mat`.

If a file starts with perXtYr, then the poses' translations have been perturbed with (X/10)% zero-mean Gaussian noise and their rotations have been perturbed (Y/10)% zero-mean Gaussian noise.

## Data in File
The results are stored in a python `dict()`. The keys for the relevent data are:

- **dual_time**: The time needed to solve the dual problem
- **dual_gt_value**: The cost function value using the ground truth solution
- **dual_primal**: The primal cost function value using the dual optimal primal solution
- **dual_gap**: The gap between the dual cost function value and the primal cost function value
- **dual_trans_error**: The magnitude of the translation error between the dual estimated and ground truth extrinsic calibration
- **dual_rot_error**: The magnitude of the rotation error between the dual estimated and ground truth extrinsic calibration
- **dual_alpha_error**: The absolute difference between the dual estimated and ground truth scale
- **rel_time**: The time needed to solve the SDP relaxation problem
- **rel_gt_value**: The cost function value using the ground truth solution
- **rel_primal**: The primal cost function value using the SDP relaxation optimal primal solution
- **rel_gap**: The gap between the SDP relaxation cost function value and the primal cost function value
- **rel_trans_error**: The magnitude of the translation error between the SDP relaxation estimated and ground truth extrinsic calibration
- **rel_rot_error**: The magnitude of the rotation error between the SDP relaxation estimated and ground truth extrinsic calibration
- **rel_alpha_error**: The absolute difference between the SDP relaxation estimated and ground truth scale.

### Important Notes
- If N runs were performed, then these keys will return arrays of length N. Otherwise, they will return a scalar value
- If the solver did not find a solution, then all of its associated values in the results are set to 1000
