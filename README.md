# MyoMesh
## Pre-Requisites
- FEniCS 2019.1.0.
- Gmsh
- LDRB (conda)
- meshio
- h5py 
- Scipy
- CMake
- VTK (libvtk7-dev)
- [hexa-mesh-from-VTK](https://github.com/rsachetto/hexa-mesh-from-VTK.git): This repository is necessary for the generation of hexahedral meshes from VTK files. It will be cloned during the Configuration.
  

## Configuration
  ```sh
    bash config.sh
  ```

## Description parameters
- `-i`: Path to the file with heart meshes.
- `-o`: Name for the output file.
- `-dx`, `-dy`, and `-dz`: refer to the discretization for the `.vtu`. Conventionally, we use the value of 0.5.
- `--alpha_endo_lv`: Fiber angle on the left ventricle (LV) endocardium. Default value is 30°.
- `--alpha_epi_lv`: Fiber angle on the left ventricle (LV) epicardium. Default value is -30°.
- `--beta_endo_lv`: Sheet angle on the left ventricle (LV) endocardium. Default value is 0°.
- `--beta_epi_lv`: Sheet angle on the left ventricle (LV) epicardium. Default value is 0°.
- `--alpha_endo_sept`: Fiber angle on the septum endocardium. Default value is 60°.
- `--alpha_epi_sept`: Fiber angle on the septum epicardium. Default value is -60°.
- `--beta_endo_sept`: Sheet angle on the septum endocardium. Default value is 0°.
- `--beta_epi_sept`: Sheet angle on the septum epicardium. Default value is 0°.
- `--alpha_endo_rv`: Fiber angle on the right ventricle (RV) endocardium. Default value is 80°.
- `--alpha_epi_rv`: Fiber angle on the right ventricle (RV) epicardium. Default value is -80°.
- `--beta_endo_rv`: Sheet angle on the right ventricle (RV) endocardium. Default value is 0°.
- `--beta_epi_rv`: Sheet angle on the right ventricle (RV) epicardium. Default value is 0°.

## Running
```sh
conda activate fenicsproject
```
```sh
python3 -i path_mesh -o output_file_name -dx dx -dy dy -dz dz
```
## Running example
```sh
python3 main.py -i ./patient1.msh -o output -dx 0.5 -dy 0.5 -dz 0.5
``` 

``` sh
python3 main.py -i ./patient1.msh -o output -dx 0.5 -dy 0.5 -dz 0.5 --alpha_endo_lv 30 --alpha_epi_lv -30 --beta_endo_lv 0 --beta_epi_lv 0 --alpha_endo_sept 60 --alpha_epi_sept -60 --beta_endo_sept 0 --beta_epi_sept 0 --alpha_endo_rv 80 --alpha_epi_rv -80 --beta_endo_rv 0 --beta_epi_rv 0
``` 