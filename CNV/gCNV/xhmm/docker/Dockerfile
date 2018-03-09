FROM gcc:5.1.0
WORKDIR /root
ADD GenomeAnalysisTK-3.8-0.jar /root/
ADD statgen-xhmm-cc14e528d909/ /root/statgen-xhmm-cc14e528d909 
ADD plinkseq-0.10/ /root/plinkseq-0.10 

RUN apt-get update
RUN apt-get install -y liblapack-dev

RUN cd statgen-xhmm-cc14e528d909; make

RUN echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee /etc/apt/sources.list.d/webupd8team-java.list
RUN echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee -a /etc/apt/sources.list.d/webupd8team-java.list

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886
RUN apt-get update
RUN echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections
RUN echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections

RUN apt-get install -y oracle-java8-installer
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle
ENV PATH $PATH:/usr/lib/jvm/java-8-oracle/jre/bin:/usr/lib/jvm/java-8-oracle/bin

EXPOSE 80
