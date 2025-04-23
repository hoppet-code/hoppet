The Python interface can be compiled with CMake

```
mkdir build; cd build
cmake ..
make [-j]
```
It is compiled by the main CMake file by default.


Alternatively, to generate the swig interface run

`swig -python  -c++ -I../src hoppet_v1.i`

And then the python module is generated with

`python setup.py sdist build_ext --inplace`

To install the package locally run

`pip install dist/hoppet_v1-1.3.0.tar.gz`

The hoppet interface can be loaded in python through

`import hoppet_v1`

For an example of how to use the interface, take a look at `tabulation_example.py`.
