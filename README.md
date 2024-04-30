# MoDVE

MoF3D generates the 3D forest for epiphytes to inhabit, growing across multiple timesteps.

MoDVE scripts:

A1- converts the MoF3D output into readable microhabitat matrices for each timestep for the epiphytes to inhabit. The epiphytes depend on the amount of substrate (i.e. branch) and the amount of light in each voxel of the microhabitat matrix.

A2- generates several pools of epiphyte species with random traits

A3- distributes these randomly generated epiphytes in forest in the initial timestep.

A4- simulates these communities of random species growing in the forest microhabitats in the subsequent timesteps

B1- selects for viable/realistic species to include in the final simulation- i.e. species that won't heavily dominate the community, or rapidly go extinct. It then generates new pools of these realistic species to use.

B2- distributes these selected species in the forest in the initial timestep.

B3- simulates these communities growing the forest microhabitats of the subsequent timesteps.
