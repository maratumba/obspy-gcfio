FROM centos:centos7

MAINTAINER Lion Krischer

RUN rpm -Uvh http://download.fedoraproject.org/pub/epel/7/x86_64/e/epel-release-7-5.noarch.rpm
# Can fail on occasion.
RUN yum -y upgrade || true
RUN yum install -y gcc numpy scipy python-matplotlib python-sqlalchemy python-lxml python-mock python-basemap python-basemap-data python-pip
RUN easy_install -U setuptools pip
RUN pip install flake8 suds-jurko future
