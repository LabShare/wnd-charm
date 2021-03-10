FROM ubuntu:18.04  

# Setting Language and Encoding
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y locales \
    && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

# Installing  packages for the build
RUN apt-get -y update  && apt-get install -y \
build-essential \
cmake \
libboost-all-dev \
libxalan-c-dev \
libqt5xmlpatterns5-dev \
libfmt-dev \
doxygen \
libgtest-dev \
libqt5svg5-dev \
libglm-dev \
libeigen3-dev  \
libtiff5-dev \
libpng-dev \
g++ \
git \
python3-pip \
python3.6-dev \
software-properties-common \
wget \
zip \
python3-gdbm \
make


# Installing virtual env
RUN apt-get update && apt-get install \
  -y --no-install-recommends python3 \
	python3-virtualenv 

# Setting Virtual env
ENV VIRTUAL_ENV=/opt/venv 
RUN python3 -m virtualenv --python=/usr/bin/python3 $VIRTUAL_ENV


#Resolving issue with Gtest
 RUN . /opt/venv/bin/activate && \
 cd /usr/src/gtest && \
 cmake CMakeLists.txt && \
 make  && \
 cp *.a /usr/lib 


ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tmp/ome-files/out/lib:/home/fftw/lib:/usr/lib/x86_64-linux-gnu

 # Setting WorkDir , downloading the repo ,building and installing the make files
 WORKDIR  /tmp/ 
 RUN . /opt/venv/bin/activate && \
 git clone https://github.com/pbsudharsan/ome-files.git && \ 
 cd  ome-files && \
 pip3 install -r requirements.txt && \
 git submodule update --init --recursive && \
 cmake -H. -BRelease -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./out  && \
 cmake --build Release --target install  && \
 ldconfig -v

 
WORKDIR /home/fftw
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz
RUN tar xfz fftw-3.3.8.tar.gz
RUN rm fftw-3.3.8.tar.gz
WORKDIR /home/fftw/fftw-3.3.8/
RUN ./configure --prefix=/home/fftw
RUN make install

COPY . /home/
WORKDIR /home/
ENV LDFLAGS="-L/home/fftw/lib -L/tmp/ome-files/out/lib -L/usr/lib/x86_64-linux-gnu"
ENV CPPFLAGS="-I/home/fftw/include  -I/tmp/ome-files/out/include"

RUN ./configure --prefix=/home/ LIBS='-L/tmp/ome-files/out/lib -lome-files -lome-xml -lome-xalan-util -lome-common -lome-xerces-util -L/usr/lib/x86_64-linux-gnu/ -lboost_iostreams  -lboost_system -lboost_filesystem -lxerces-c-3.2' CXXFLAGS='-g -O2 -fopenmp -I/tmp/ome-files/out/include '

RUN make install

RUN mkdir /home/input 
RUN mkdir /home/output 
RUN mkdir /home/mask 

ENTRYPOINT ["/home/bin/wndchrm"]

#sudo docker build -t labshare/polus-wnd-charm-plugin:0.2.0 .

#sudo docker run -v /home/maghrebim2/Work/WND-CHARM/ROI/MultiThreading/Working_Dir2/Data:/home/input -v /home/maghrebim2/Work/WND-CHARM/ROI/MultiThreading/Working_Dir2/LabeledData:/home/mask -v /home/maghrebim2/Work/WND-CHARM/ROI/MultiThreading/Working_Dir2/Output:/home/output  labshare/polus-wnd-charm-plugin:0.2.0 --DataPath /home/input --LabeledData /home/mask  --ImageTransformationName  Chebyshev-Fourier_1D_ColumnWise  --FeatureAlgorithmName PixelStatistics --output /home/output 
