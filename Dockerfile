FROM python:3.11

ARG SUNDIALS_VERSION=2.7.0

# System deps
RUN apt-get update && apt-get -y install cmake liblapack-dev libsuitesparse-dev libhypre-dev curl git make vim bash-completion
RUN cp -v /usr/lib/x86_64-linux-gnu/libblas.so /usr/lib/x86_64-linux-gnu/libblas_OPENMP.so

# Python packages
RUN pip install Cython numpy scipy matplotlib setuptools==69.1.0

# SuperLU MT 4.0.1
RUN cd /tmp && \
    curl -fSsL https://github.com/xiaoyeli/superlu_mt/archive/refs/tags/v4.0.1.tar.gz | tar xz && \
    cd superlu_mt-4.0.1 && \
    cmake -Denable_examples=OFF -Denable_tests=OFF -DPLAT="_OPENMP" \
          -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_INSTALL_LIBDIR=lib \
          -DSUPERLUMT_INSTALL_INCLUDEDIR=include . && \
    make -j4 && make install

# Sundials (version-parameterized; v2.7.0 needs CMakeLists patch)
RUN git clone --depth 1 -b v${SUNDIALS_VERSION} https://github.com/LLNL/sundials /tmp/sundials && \
    cd /tmp/sundials && \
    if [ "${SUNDIALS_VERSION}" = "2.7.0" ]; then \
        echo "target_link_libraries(sundials_idas_shared lapack blas superlu_mt_OPENMP)" >> src/idas/CMakeLists.txt; \
        echo "target_link_libraries(sundials_kinsol_shared lapack blas superlu_mt_OPENMP)" >> src/kinsol/CMakeLists.txt; \
    fi && \
    mkdir build && cd build && \
    cmake -LAH \
        -DSUPERLUMT_BLAS_LIBRARIES=blas \
        -DSUPERLUMT_LIBRARIES=blas \
        -DSUPERLUMT_INCLUDE_DIR=/usr/include \
        -DSUPERLUMT_LIBRARY=/usr/lib/libsuperlu_mt_OPENMP.a \
        -DSUPERLUMT_THREAD_TYPE=OpenMP \
        -DCMAKE_INSTALL_PREFIX=/usr \
        -DSUPERLUMT_ENABLE=ON \
        -DLAPACK_ENABLE=ON \
        -DEXAMPLES_ENABLE=OFF \
        -DEXAMPLES_ENABLE_C=OFF \
        -DBUILD_STATIC_LIBS=OFF \
        -DSUNDIALS_INDEX_SIZE=32 \
        .. && \
    make -j4 && make install

ARG PYTHON_VENV=/src/.venv
ENV PATH=${PYTHON_VENV}/bin:$PATH
WORKDIR /src
