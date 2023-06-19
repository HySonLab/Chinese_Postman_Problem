# Implementation in C++ for the load-dependent Chinese postman problem
## Author: Dr. Truong Son Hy
## Copyright 2023

### Data generation (```data_generation/```)

* Data generation is based on paper "The Chinese Postman Problem with Load-Dependent Costs" https://pubsonline.informs.org/doi/abs/10.1287/trsc.2017.0774
* Implementation to generate Eulerian cases: ```gen_eulerian.cpp```.
* Implementation to generate Rural Postmam Problem (RPP) cases: ```gen_rural.cpp```.
* In the ```data/``` directory: Prefix `E` denotes the Eulerian cases. Prefix `C` denotes the cases based on Christofides et al. Prefix `H` denotes the cases based on Hertz et al.

