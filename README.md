Personal code for principal component analysis and diffusion map
[examples](#examples).  Specifically made to test the idea on some well-known
data sets, but it wouldn't take much to modify the source for use with whatever
data set or distance metric you desire.

## Compilation

    $ make

A library is compiled with the classes needed for the main program and the main
program links to that. The main program requires
[json-fortran](https://github.com/jacobwilliams/json-fortran). LAPACK is
required for the library to calculate the eigenvectors and eigenvalues of
various matrices.

## Running

Modify `dmap.json`. Then do:

    $ ./run dmap.json

You can also run principal component analysis using the following file:

    $ ./run pca.json

`bandwidth.json` is for running the program iteratively over different bandwidth
values. See Figure S1 in [this
document](https://www.pnas.org/cgi/doi/10.1073/pnas.1003293107) for what I was
going for with this. This would more helpful for analyzing simulation data, but
the main program is not set up for that.

## Extras

The `extras` folder contains the source code of two programs to aid in
generating example data sets. No configuration files are provided, so you will
need to edit the source.

## Examples

Compare the following results with those found [in this
paper](http://wireilla.com/papers/ijfcst/V4N6/4614ijfcst06.pdf), specifically in
Section 3.1. Note that my value of `bandwidth` is the square of what they
call `sigma` (I am not squaring the denominator of the Gaussian kernel in my
code).

### Cluster of points

Colors indicate where points are in relationship to axis with greatest variance.

#### Original data
![Cluster dataset](examples/cluster/cluster.png)

#### Principal component analysis
![PCA](examples/cluster/pca.png)

#### Diffusion maps
![Diffusion map (eps = 1)](examples/cluster/dmap-1.png)
![Diffusion map (eps = 10)](examples/cluster/dmap-10.png)
![Diffusion map (eps = 100)](examples/cluster/dmap-100.png)

### Swiss roll

Colors indicate where points are in relationship to the center of the swiss
roll.

#### Original data
![Swiss roll dataset](examples/swiss-roll/swissroll.png)

#### Principal component analysis
![PCA](examples/swiss-roll/pca.png)

#### Diffusion maps
![Diffusion map (eps = 1)](examples/swiss-roll/dmap-1.png)
![Diffusion map (eps = 10)](examples/swiss-roll/dmap-10.png)
![Diffusion map (eps = 100)](examples/swiss-roll/dmap-100.png)

### Punctured sphere

Colors indicate where points are in relationship to axis that goes through the
holes in the sphere.

#### Original data
![Punctured sphere dataset](examples/punctured-sphere/puncsphere.png)

#### Principal component analysis
![PCA](examples/punctured-sphere/pca.png)

#### Diffusion maps
![Diffusion map (eps = 1)](examples/punctured-sphere/dmap-1.png)
![Diffusion map (eps = 10)](examples/punctured-sphere/dmap-10.png)
![Diffusion map (eps = 100)](examples/punctured-sphere/dmap-100.png)

