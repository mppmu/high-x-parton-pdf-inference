FROM fedora:39
RUN   dnf -y update && dnf -y install dnf5  && dnf -y clean all 
RUN   dnf5 -y install julia python gcc-c++ gcc-gfortran && dnf5 -y clean all
