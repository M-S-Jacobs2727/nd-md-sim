# n-Dimensional Molecular Dynamics Simulator

My plans for this moving forward (effective 2023-12-24) are to start development again, focusing on
a single-node version using shared memory parallelism with one or two potentials and no molecular
topology. Then I will add, in some order, GPU acceleration (single card for shared-memory), more
complicated potentials, and bond, angle, dihedral, and improper potentials, as well as the possibility
of writing custom potentials very easily. 

Note: move Atoms struct into separate file, split main file apart, and make everything into functions.