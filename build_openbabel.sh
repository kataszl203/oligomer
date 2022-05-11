#!/bin/sh
PYTHON_EXECUTABLE=$(which python3)
echo $PYTHON_EXECUTABLE
SCRIPT=$(readlink -f "$0")
WORKDIR=$(dirname "$SCRIPT")
mkdir -p $WORKDIR/external && cd $WORKDIR/external

if [ ! -d "./openbabel" ]; then
  git clone https://github.com/openbabel/openbabel.git
fi

cd openbabel &&\
  git checkout f3ed2a9a5166dbd3b9ce386e636a176074a6c34c &&\
  git reset --hard HEAD &&\
  git apply $WORKDIR/openbabel/openbabel-fix.patch &&\
  mkdir -p build &&\ 
  cd build

cmake ..\
  -DPYTHON_EXECUTABLE=$PYTHON_EXECUTABLE\
  -DPYTHON_BINDINGS=ON\
  -DRUN_SWIG=ON\
  -DWITH_MAEPARSER=off\
  -DCMAKE_INSTALL_PREFIX=$WORKDIR/external/build &&\
  nproc=$(getconf _NPROCESSORS_ONLN) &&\
  make -j $(( nproc > 2 ? nproc - 2 : 1 )) &&\
  make install

rm -rf $WORKDIR/external/build/lib/python
mv $WORKDIR/external/build/lib/python* $WORKDIR/external/build/lib/python

echo
echo "build success!"
echo
