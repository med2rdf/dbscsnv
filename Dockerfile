FROM ruby:2.5

RUN gem install rdf-vocab
RUN gem install rdf-turtle

RUN mkdir /work
WORKDIR /work

ADD mkRDFdbscSNV.rb dbscsnv.conf /work/

ENV RUBYOPT -EUTF-8

