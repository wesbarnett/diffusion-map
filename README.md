Personal code for diffusion map examples. Only suitable for 3-dimensional data
at this point, so its really just useful for toy data sets. However, it
wouldn't take much to modify the source for use with whatever data set or
distance metric you desire.

## Examples

Compare the following results with those found [in this
paper](http://wireilla.com/papers/ijfcst/V4N6/4614ijfcst06.pdf), specifically in
Section 3.1. Note that my value of `bandwidth` is what the square of what they
call `sigma` (I am not squaring the denominator of the Gaussian kernel in my
code).

### Swiss roll

![Swiss roll dataset](examples/swiss-roll/swissroll.png)
![Diffusion map (eps = 1)](examples/swiss-roll/dmap-1.png)
![Diffusion map (eps = 10)](examples/swiss-roll/dmap-10.png)
![Diffusion map (eps = 100)](examples/swiss-roll/dmap-100.png)

### Punctured sphere

![Punctured sphere dataset](examples/punctured-sphere/puncsphere.png)
![Diffusion map (eps = 1)](examples/punctured-sphere/dmap-1.png)
![Diffusion map (eps = 10)](examples/punctured-sphere/dmap-10.png)
![Diffusion map (eps = 100)](examples/punctured-sphere/dmap-100.png)

## Compilation

    $ make

## Running

Modify `run.json`. Then do:

    $ ./run run.json

`bandwidth.json` is for running the program iteratively over different bandwidth
values. See Figure S1 in [this
document](https://www.pnas.org/cgi/doi/10.1073/pnas.1003293107) for what I was
going for with this.

## Extras

The `extras` folder contains the source code of two programs to aid in
generating example data sets. No configuration files are provided, so you will
need to edit the source.

