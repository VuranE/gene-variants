# gene-variants

# Discovering Gene Variants From Sequencing Data

This is a student project developed for the Bioinformatics 1 course at the Faculty of Electrical Engineering and Computing (FER), University of Zagreb. Course information can be found here: https://www.fer.unizg.hr/en/course/enbio1.

The goal of the project is to identify gene variants from sequencing reads using clustering techniques based on sequence similarity. The SPOA (Simd Partial Order Alignment) library is used to compute consensus sequences for each cluster, representing distinct gene variants.

How to install:

<pre>
git clone https://github.com/VuranE/gene-variants.git
cd gene-variants 
mkdir -p external
cd external
git clone https://github.com/rvaser/spoa.git
cd ..
sudo apt update
sudo apt install build-essential cmake
mkdir -p build
cd build
cmake ..
make
</pre>

Running the program:
 `./project <path-to-sample-files-directory> <path-to-ground-truth-files>`
