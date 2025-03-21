# Magnetically Actuated Capsule Control – Reproduction of Pietro Valdastri's Work

This project is a reproduction of the work conducted by **Pietro Valdastri's research group**. It simulates the control of a robotic arm to drive a capsule using magnetic fields. The simulation environment enables users to experiment with **magnetic capsule manipulation** in a controlled setting.

## Installation Requirements

To run this simulation, the following software is required:

1. **MATLAB** (with the **Robotics System Toolbox**)  
   - Install the Robotics System Toolbox from:  
     [MathWorks Robotics System Toolbox](https://www.mathworks.com/products/robotics.html)
2. **CoppeliaSim** (formerly known as V-REP)  
   - Used for the simulation of the robotic arm and magnetic interaction.
  
## Usage

1. Open **CoppeliaSim** and run `akukamag_ipm.ttt`.  
2. In **MATLAB**, run `EndoCtl_ver1.m` to start the control process.  

## Reference

This MATLAB implementation is based on the following research paper:  

@article{taddese2018enhanced,  
title={Enhanced real-time pose estimation for closed-loop robotic manipulation of magnetically actuated capsule endoscopes},  
author={Taddese, Addisu Z and Slawinski, Piotr R and Pirotta, Marco and De Momi, Elena and Obstein, Keith L and Valdastri, Pietro},  
journal={The International Journal of Robotics Research},  
volume={37},  
number={8},  
pages={890--911},  
year={2018},  
publisher={SAGE Publications Sage UK: London, England}  
}

## Further Reading  

For research on **manipulator control under floating-base disturbances**, consider:  

@article{xu2024confidence,  
  title={Confidence-Aware Object Capture for a Manipulator Subject to Floating-Base Disturbances},  
  author={Xu, Ruoyu and Jiang, Zixing and Liu, Beibei and Wang, Yuquan and Qian, Huihuan},  
  journal={IEEE Transactions on Robotics},  
  year={2024},  
  publisher={IEEE}  
}

Although unrelated to this project, this paper explores techniques for robust robotic manipulation under uncertainty, which may be of interest to researchers working on robotic control.  



