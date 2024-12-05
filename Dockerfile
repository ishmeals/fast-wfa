FROM ubuntu:20.04
ENV DEBIAN_FRONTEND="noninteractive" TZ="America/New_York"
RUN apt-get update 
RUN apt install -y lsb-release wget software-properties-common gnupg

RUN wget https://apt.llvm.org/llvm.sh
RUN chmod +x llvm.sh
RUN ./llvm.sh 19 all

RUN apt-get install -y git build-essential pip
#RUN apt-get install -y curl zip unzip tar pkg-config linux-libc-dev

RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install cmake
RUN python3 -m pip install ninja 


WORKDIR ~/
RUN mkdir fast-wfa
RUN mkdir fast-wfa/src
RUN mkdir fast-wfa/include
RUN mkdir fast-wfa/tests
#RUN mkdir SynthesisSearch/vcpkg_overlays
RUN mkdir fast-wfa/.git
COPY CMakeLists.txt CMakePresets.json fast-wfa
COPY src fast-wfa/src
COPY include fast-wfa/include
COPY tests fast-wfa/tests
#COPY vcpkg_overlays SynthesisSearch/vcpkg_overlays
COPY .git SynthesisSearch/.git
WORKDIR fast-wfa
#RUN git submodule init
#RUN git submodule update
#RUN update-alternatives --install /usr/bin/c++ c++ /usr/bin/clang++ 50

RUN cmake --preset=linux-release
WORKDIR ./out/build/linux-release
RUN ninja
WORKDIR ./bin

CMD ["./wfa_tool"]