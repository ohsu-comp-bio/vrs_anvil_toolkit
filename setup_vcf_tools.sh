cd ~
curl -LJO https://github.com/vcftools/vcftools/tarball/master

VCF_TOOLS_TAR=$(ls -1t vcftools*.tar.gz | head -n 1)
tar -xzvf $VCF_TOOLS_TAR
rm $VCF_TOOLS_TAR

VCF_TOOLS_DIR=$(ls -1td vcftools* | head -n 1)
cd $VCF_TOOLS_DIR

./autogen.sh
./configure
make
# TODO: test with sudo
sudo make install

# TODO: test executable exists in the dir or whatever