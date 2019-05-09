FROM alpine:3.6
MAINTAINER Laurence Yang <lyang@eng.ucsd.edu>

# TODO: Add the source code
# ADD src.tar.gz /source/

RUN echo "@community http://dl-cdn.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories

ENV PACKAGES="\
	dumb-init \
	musl \
	linux-headers \
	build-base \
	gfortran \
	bash \
	git \
	ca-certificates \
	openblas-dev@community \
	python2 \
	python2-dev \
	py-setuptools \
"

RUN echo \
	# Add packages
	&& apk add --no-cache $PACKAGES \
	# symlinks
	&& if [[ ! -e /usr/bin/python ]]; 			then ln -sf /usr/bin/python2.7 /usr/bin/python; fi \
	&& if [[ ! -e /usr/bin/python-config ]]; 	then ln -sf /usr/bin/python2.7-config /usr/bin/python-config; fi \
	&& if [[ ! -e /usr/bin/easy_install ]]; 	then ln -sf /usr/bin/easy_install-2.7 /usr/bin/easy_install; fi \
	# Install and upgrade pip
	&& easy_install pip \
	&& pip install --upgrade pip \
	&& if [[ ! -e /usr/bin/pip ]]; then ln -sf /usr/bin/pip2.7 /usr/bin/pip; fi \
	&& pip install scipy \
	&& pip install pandas \
	# TODO: copy the source code
	&& make clean \
	&& make \
	&& cd /source \
	&& mkdir python_modules \
	# TODO: copy the compiled binaries
	&& python setup.py develop --user \
	&& echo

# mount volume
ENTRYPOINT ["/usr/bin/dumb-init"]
