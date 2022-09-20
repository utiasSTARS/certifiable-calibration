# certifiable-calibration
Certifiably globally optimal extrinsic calibration for sensors providing egomotion estimates. 

<img src="https://raw.githubusercontent.com/utiasSTARS/certifiable-calibration/master/calibration_high_level.png" width="500px"/>


## MATLAB Code for Egomotion Sensor Calibration
See the [`matlab/`](https://github.com/utiasSTARS/certifiable-calibration/tree/master/matlab) directory for instructions on how to run the algorithm from [our 2019 IEEE RA-L publication](https://arxiv.org/pdf/1809.03554.pdf). The [`data/`](https://github.com/utiasSTARS/certifiable-calibration/tree/master/data) directory contains files with experimental results from our paper.


## Python Code for Monocular Hand-Eye Calibration
See the [`python/extrinsic_calibration`](https://github.com/utiasSTARS/certifiable-calibration/tree/master/python/extrinsic_calibration) directory for instructions on how to run the algorithm from [our IEEE MFI 2020 publication](https://arxiv.org/abs/2005.08298).

##  Citation
If you use this work in your research, please cite the following papers:

```
@article{2019_Giamou_Certifiably,
  doi = {10.1109/LRA.2018.2890444},
  journal = {{IEEE} Robotics and Automation Letters},
  month = {April},
  number = {2},
  pages = {367--374},
  title = {Certifiably Globally Optimal Extrinsic Calibration from Per-Sensor Egomotion},
  url = {https://arxiv.org/abs/1809.03554},
  volume = {4},
  year = {2019}
}

@inproceedings{2020_Wise_Certifiably,
  address = {Karlsruhe, Germany},
  author = {Emmett Wise and Matthew Giamou and Soroush Khoubyarian and Abhinav Grover and Jonathan Kelly},
  booktitle = {Proceedings of the {IEEE} International Conference on Multisensor Fusion and Integration {(MFI)}},
	doi = {10.1109/MFI49285.2020.9235219},
	title = {Certifiably Optimal Monocular Hand-Eye Calibration},
	url = {https://arxiv.org/abs/2005.08298},
	year = {2020}
}
```
