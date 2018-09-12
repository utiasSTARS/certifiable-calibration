# certifiable-calibration
Certifiably globally optimal extrinsic calibration for sensors providing egomotion estimates. 

<img src="https://raw.githubusercontent.com/utiasSTARS/certifiable-calibration/master/calibration_high_level.png" width="500px"/>


## Installation and Dependencies 

All code and experiments were developed in MATLAB R2017b.

We use the [CVX](http://cvxr.com/cvx/) modelling language, which is free for academic usage, for all convex optimization. The default SDPT3 solver that ships with CVX was used.

Our plotting script makes use of the [`subaxis`](https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot) function to create subplots. 

Our experiments made use of the code and dataset available [here](http://jbrookshire.com/projects_3dcalib.htm), which is not freely available but can be requested from the author for academic purposes. Much of the data loaded and used in our results comes from applications of that code.

## Usage 
All code is found in the folder `matlab/`. Be sure to add this folder and its subfolders to your MATLAB path.

### Plotting Script
The code for plots appearing in our [paper](https://arxiv.org/pdf/1809.03554.pdf) listed below can be run with the script `matlab/plot_ral_figures.m`. This script uses the `.mat` files containing experimental results in `data/`. The data files are addressed relatively, so be sure to run the script from within the `matlab/` folder.

### Calibration Example Script
The script `matlab/example_egomotion_calibration.m` shows a sample usage of the function `matlab/egomotion_calibration/egomotion_calibration.m` on simulated noisy data. 

<img src="https://raw.githubusercontent.com/utiasSTARS/certifiable-calibration/master/sensor_egomotion.png" width="400px"/>


## Citation
If you use any of this code in your work, please cite the [relevant publication](https://arxiv.org/pdf/1809.03554.pdf): 

```bibtex
@article{giamou2017talk,
  title   = {Certifiably Globally Optimal Extrinsic Calibration from Per-Sensor Egomotion},
  author  = {Giamou, Matthew and Ma, Ziye and Peretroukhin, Valentin and Kelly, Jonathan},
  journal = {{IEEE} Robotics and Automation Letters (submitted)}
  year    = {2018}
}
```
