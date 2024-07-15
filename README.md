# ARapidTb
Targeted sputum sequencing for rapid and comprehensive drug resistance of Mycobacterium tuberculosis

**To run with Docker**

``git clone https://github.com/jade-nhri/ARapidTb.git``

``cd ARapidTb``

``docker build -t "rapidtb:v1" ./``

``docker run --runtime=nvidia -h rapidtb --name rapidtb -i -t -v /:/MyData rapidtb:v1 /bin/bash``


Installation
------------
**Installation from source**

``cd /opt``

``git clone https://github.com/jade-nhri/ARapidTb.git``

``cd ARapidTb``

``chmod +x *.py``

``export PATH="$PATH:/opt/ARapidTb"``


## Dependencies

- [pyspoa-0.0.3](https://github.com/nanoporetech/pyspoa)
- [medaka-1.5.0](https://github.com/nanoporetech/medaka)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools-1.13](http://github.com/samtools/)
- [bcftools](https://github.com/samtools/bcftools)
- [seqkit-v2.8.0](https://github.com/shenwei356/seqkit)

Quick start
------------
funfilter.py inputdir /opt/ARapidTb/Tbrefs.fasta outputdir
preduct.py indir /opt/ARrapidTb/Tbrefs.fasta outputdir
