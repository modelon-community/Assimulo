.PHONY: build build-dev-image test shell

SUNDIALS_VERSION ?= 2.7.0
DOCKER_IMAGE := assimulo-dev-sundials_$(SUNDIALS_VERSION)

IN_DOCKER_IMG := $(shell test -f /.dockerenv && echo 1 || echo 0)

FORTRAN_FLAGS := "-std=legacy"
SETUPTOOLS_JFLAG=-j$(shell nproc)

define _run
@if [ $(IN_DOCKER_IMG) -eq 1 ]; then \
$(1);\
else \
docker run \
--rm $(2) \
-v $(CURDIR):/src \
$(DOCKER_IMAGE) \
$(1); \
fi
endef

build-dev-image:
	docker build --build-arg SUNDIALS_VERSION=$(SUNDIALS_VERSION) -t $(DOCKER_IMAGE) .

.venv:
	$(call _run, python3.11 -m venv .venv --system-site-packages)
	$(call _run, pip install pytest)

build: .venv
	$(call _run, python3.11 setup.py build_ext ${SETUPTOOLS_JFLAG} install --sundials-home=/usr --blas-home=/usr/lib/x86_64-linux-gnu/ --lapack-home=/usr/lib/x86_64-linux-gnu/ --superlu-home=/usr --extra-fortran-compile-flags=$(FORTRAN_FLAGS))

test: build
	$(call _run, pytest)

shell:
	$(call _run, /bin/bash,-it)
