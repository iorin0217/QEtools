FROM ubuntu:18.04

# sh doesn't support 'source' command
SHELL ["/bin/bash", "-c"]

RUN apt update &&\
    apt -y upgrade &&\
    DEBIAN_FRONTEND=noninteractive apt -y install wget environment-modules gnupg rsh-client tclsh build-essential

# install pgi
# --with-cma="no" or sudo echo 0 > /proc/sys/kernel/yama/ptrace_scope will be applied
ENV PGI_VERSION 19.4
ENV PGI_RELEASE 2019
ADD pgilinux-2019-194-x86-64.tar.gz /tmp
RUN PGI_SILENT=true\
    PGI_ACCEPT_EULA=accept\
    PGI_INSTALL_NVIDIA=false\
    PGI_INSTALL_MANAGED=true\
    PGI_INSTALL_AMD=false\
    PGI_INSTALL_JAVA=false\
    PGI_INSTALL_MPI=true\
    PGI_MPI_GPU_SUPPORT=false\
    /tmp/install &&\
    echo "source /usr/share/modules/init/bash" >> /etc/bash.bashrc &&\
    echo "module load pgi openmpi" >> /etc/bash.bashrc &&\
    echo "/opt/pgi/modulefiles" >> /etc/environment-modules/modulespath &&\
    rm -rf /tmp/*

# install intel-mkl and intel-mpi
ENV MKL_VERSION 2019.3-062
ADD https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB /tmp
RUN apt-key add /tmp/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB &&\
    rm /tmp/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB &&\
    echo "deb https://apt.repos.intel.com/mkl all main" > /etc/apt/sources.list.d/intel-mkl.list &&\
    apt update &&\
    apt install -y intel-mkl-${MKL_VERSION} &&\
    echo /opt/intel/mkl/lib/intel64 > /etc/ld.so.conf.d/intel_mkl.conf &&\
    ldconfig

# build intel mkl fftw with non-llvm pgi
# ( build with llvm pgi fails)
RUN cd /opt/intel/mkl/interfaces/fftw3xc &&\
    PATH=/opt/pgi/linux86-64-nollvm/${PGI_VERSION}/bin:${PATH} make compiler=pgi libintel64 &&\
    mv libfftw3xc_pgi.a /opt/intel/mkl/lib/intel64

# build quantum-espresso with intel-mkl
ENV QE_VERSION 6.4.1
ADD https://gitlab.com/QEF/q-e/-/archive/qe-${QE_VERSION}/q-e-qe-${QE_VERSION}.tar.gz /tmp
RUN cd /tmp &&\
    tar xf q-e-qe-${QE_VERSION}.tar.gz &&\
    cd /tmp/q-e-qe-${QE_VERSION} &&\
    source /usr/share/modules/init/bash &&\
    module load pgi openmpi &&\
    source /opt/intel/bin/compilervars.sh intel64 &&\
    ./configure CPP=cpp CC=pgcc F90=pgf90 MPIF90=mpifort\
    BLAS_LIBS="-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"\
    LAPACK_LIBS="-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"\
    SCALAPACK_="-L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64"\
    FFT_LIBS="/opt/intel/mkl/lib/intel64/libfftw3xc_pgi.a" &&\
    make pwall &&\
    make install

RUN apt-get update -y && \
    apt-get install -y vim tmux git ssh &&\
    apt-get install -y python3 python3-pip &&\
    pip3 install numpy scipy pandas matplotlib plotly pymatgen seekpath ase

# add account
RUN groupadd -g 2019 Keisan
RUN useradd -u 2019 -g 2019 -d /home/Keisan --create-home --shell /usr/bin/zsh Keisan && \
   echo "Keisan ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
   chown -R Keisan:Keisan /home/Keisan
ENV HOME /home/Keisan

RUN rm -rf /var/lib/apt/lists/* &&\
    rm /tmp/* -rf