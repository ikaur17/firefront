# README #

As a first step, execute the following command:
```
#!bash
$ chmod +x ./move ./clean ./script_*.sh
```
then, in order to compile the program,

* with OpenMP support
```
#!bash
$ make clean
$ make
```
* without OpenMP support:
```
#!bash
$ make clean
$ make flag_omp=no
```

### To run the first test case ###

Execute the following commands in a terminal:
```
#!bash
$ ./script_1.sh
```
To post-process the results:

```
#!bash
$ python makefig_1.py
```

This will create and place in the current directory the following files: 

* fig_a.pdf, fig_b.pdf, fig_c.pdf, fig_d.pdf, fig_e.pdf, fig_f.pdf
