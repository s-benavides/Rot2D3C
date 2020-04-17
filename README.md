# Two-dimensional, three-component turbulence including rotation (Rot2D3C)

This is a two-dimensional pseudospectral code written by Pablo Mininni and modified by Santiago J. Benavides to include the out-of-plane velocity, which is coupled to the vertical vorticity with a rotation in the y-direction.

The system of equations being solved is the following:

<img src="https://render.githubusercontent.com/render/math?math=\partial_t \omega %2B [\omega,\psi] = 2\Omega \partial_y v_z %2B \nu \nabla^2 \omega - \nu^{%2B} \nabla^{-4} \omega %2B f_\omega">

<img src="https://render.githubusercontent.com/render/math?math=\partial_t v_z %2B [v_z,\psi] = 2\Omega \partial_y \psi %2B \nu \nabla^2 v_z - \nu^{%2B} \nabla^{-4} v_z %2B f_{v_z}">

where omega is the vorticity in the out-of-plane direction, psi is the corresponding streamfunction, and vz is the out-of-plane velocity.

