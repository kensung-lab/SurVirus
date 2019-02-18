unzip htslib-1.9.zip
cd htslib-1.9
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
