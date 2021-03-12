FROM ubuntu:18.04  

# Installing  packages for the build
RUN apt-get -y update  && apt-get install -y \
build-essential \
libboost-all-dev \
libtiff5-dev \
wget \
make

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fftw/lib:/usr/lib/x86_64-linux-gnu
 
WORKDIR /home/fftw
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz
RUN tar xfz fftw-3.3.8.tar.gz
RUN rm fftw-3.3.8.tar.gz
WORKDIR /home/fftw/fftw-3.3.8/
RUN ./configure --prefix=/home/fftw
RUN make install

COPY . /home/
WORKDIR /home/
ENV LDFLAGS="-L/home/fftw/lib -L/usr/lib/x86_64-linux-gnu"
ENV CPPFLAGS="-I/home/fftw/include"
RUN ./configure --prefix=/home/ CXXFLAGS='-g -O2 -fopenmp'
RUN make install

RUN mkdir /home/input 
RUN mkdir /home/output 

ENTRYPOINT ["/home/bin/wndchrm"]

#sudo docker build -t labshare/polus-wnd-charm-plugin:0.3.0 .

#sudo docker run -v /master/Data:/home/input -v /master/Output:/home/output  labshare/polus-wnd-charm-plugin:0.3.0 train -l /home/input /home/output/output.fit 
# -l in above is optional and will force WND-CHARM to compute the long list of features (~3000)
