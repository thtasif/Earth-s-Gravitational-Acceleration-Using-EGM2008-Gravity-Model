# Earth-s-Gravitational-Acceleration-Using-EGM2008-Gravity-Model

In this project, a functon for MATLAB has been developed to compute Earth's gravitational acceleration using EGM-2008 gravity model. The function works 8 to 10 times
faster comparing to ‘gravitysphericalharmonic’ function of Aerospace toolbox of MATLAB when tested in MATLAB R2019a and using 120 degree and order of EGM-2008.
Another issue of using ‘gravitysphericalharmonic’ function is that it always takes degree and order as a single input. Hence, it is not possible to get results for
120 degree and 100 order. But using this new functon it is possible to use different values for degree and order.
