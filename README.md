# RDF tool for dbscSNV

---

## How to make RDFs of dbscSNV

`# cd /hoge`

`# mkdir data`

`# cd data`

Download data from 
https://sites.google.com/site/jpopgen/dbNSFP

Get the source code of this tool from github.

`# cd /hoge`

`# git clone https://github.com/med2rdf/dbscsnv.git`

`# cd dbscsnv`

Run docker commands.

`# docker build -t dbscsnv .`

For proxy

`#docker build -t dbscsnv --build-arg HTTPS_PROXY=http://hoge:8080 .`

Now, let's make RDFs!

`# docker run -v /hoge/data:/data dbnsfp ruby /work/mkRDFdbscSNV.rb /work/dbscsnv.conf /data/dbscSNVxxx.chrx`

- Replace dbscSNVxxx.chrx with the name of the downloaded file when you execute this command.

It outputs /hoge/data/dbscSNVxxx.chrx.ttl .
- You can replace dbscsnv.conf for customizing the rdf.

---
## Citation

This schema was drawed by
https://www.kanzaki.com/works/2009/pub/graph-draw

