unzip htslib-1.9.zip
cd htslib-1.9
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
cd ..

unzip sparsehash-sparsehash-2.0.3.zip
cd sparsehash-sparsehash-2.0.3
./configure --prefix=`pwd`
make
make install
