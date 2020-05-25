An implementation of Andreff's algorithm.

unscaled_pose_data.mat is a dataset with 0 noise. Currently the code works perfectly with this dataset.

For input_python.mat, which is a dataset accompanied with noise, the code shows a significant amount of deviation from the ground truth.
This is because the rotation matrix has not been projected onto SO3 yet in this project.