FROM ubuntu:16.04
# install
RUN apt-get update
RUN apt-get install -y apt-utils g++-4.7
RUN ln -s /usr/bin/g++-4.7 /usr/bin/g++
RUN apt-get install -y make gcc zlib1g-dev git
# cleanup
RUN 	apt-get  -y autoremove
RUN 	apt-get  -y clean
# clone
RUN git clone https://github.com/aquaskyline/16GT.git
# compile
WORKDIR 16GT
RUN make
RUN git clone https://github.com/aquaskyline/SOAP3-dp.git
WORKDIR SOAP3-dp
RUN make SOAP3-Builder
RUN make BGS-Build


