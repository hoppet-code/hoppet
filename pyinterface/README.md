To generate the swig interface run

swig -python  -c++ -I../src hoppet_v1.i

And then the python module is generated with

python setup.py sdist build_ext --inplace

To install the package locally run

pip install dist/hoppet_v1-1.3.0.tar.gz 