# CPPM: Chi-squared Progressive Photon Mapping

![result](https://bactlink.github.io/CPPM/papers_268s3.jpg)

## Introduction
This is the source code of ***CPPM: Chi-squared Progressive Photon Mapping***.
The paper can be found on the [homepage](https://bactlink.github.io/CPPM) of this project.

## Documentation
This project is based on [Mitsuba 0.6.0](https://github.com/mitsuba-renderer/mitsuba).

Our codes locate in the folder `mitsuba-CPPM/src/integrators/cppm`:
- CPPM: cppm1.cpp
- CPPM-LF: cppm2.cpp
- CPPM-GOF: cppm3.cpp

We did not modify the original codes, so feel free to copy these additional files to your own Mitsuba Renderer:
```
include/mitsuba/render/cppm*.h
src/integrators/cppm/*
src/librender/cppm*.cpp
```
Remember to modify ```src/integrators/SConscript``` to include the integrators.

## Scenes
We sort out some scenes in Mitsuba's format, which can be found [here](https://github.com/bacTlink/mitsuba-CPPM-scenes).

## BibTeX
```
@article{lin2020cppm,
  title={CPPM: Chi-squared Progressive Photon Mapping},
  author={Zehui Lin and Sheng Li and Xinlu Zeng and Congyi Zhang and Jinzhu Jia and Guoping Wang and Dinesh Manocha},
  journal={ACM Transactions on Graphics (TOG)},
  volume={39},
  number={6},
  article={240},
  year={2020},
  publisher={ACM},
  DOI = {10.1145/3414685.3417822}
}
```
