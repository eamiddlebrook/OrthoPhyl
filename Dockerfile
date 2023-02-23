FROM ubuntu:20.04
LABEL AUTHOR=earlm@lanl.gov

RUN apt-get update \
    && apt-get install -y build-essential git \
    && apt-get install -y python3 \
    && apt-get install -y wget zip bc \
    && cd /opt \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && bash miniconda.sh -bfp /opt/miniconda3 \
    && rm miniconda.sh \
    && export PATH="/opt/miniconda3/bin:$PATH" \
    && echo "PATH=\"/opt/miniconda3/bin:$PATH\"" >> /etc/environment \
    && . /opt/miniconda3/etc/profile.d/conda.sh \
    && conda activate base \
    && conda install -y -c bioconda -c conda-forge prodigal orthofinder bbmap fasttree hmmer pal2nal prodigal ete3 raxml trimal parallel r-essentials \
    && ln -s /opt/miniconda3/bin/conda /usr/bin/conda \
    && chmod --recursive a+rw /opt/miniconda3 \
    && Path_to_gits=/opt/Path_to_gits \
    && echo "Path_to_gits=/opt/Path_to_gits" >> /etc/environment \
    && mkdir $Path_to_gits \
    && cd $Path_to_gits/ \
    && git clone https://github.com/smirarab/ASTRAL.git \
    && cd ASTRAL \
    && unzip Astral.5.7.8.zip \
    && cd $Path_to_gits/ \
    && git clone https://github.com/nylander/catfasta2phyml.git \
    && git clone https://github.com/dportik/Alignment_Assessment.git \
    && cd Alignment_Assessment/ \
    && 2to3 -w Alignment_Assessment_v2.py \
    && cd /opt \
    && wget https://github.com/ParBLiSS/FastANI/releases/download/v1.33/fastANI-Linux64-v1.33.zip \
    && unzip fastANI-Linux64-v1.33.zip \
    && mkdir /opt/apps/ \
    && mv fastANI /opt/apps/ \
    && rm fastANI-Linux64-v1.33.zip \
    && cd $Path_to_gits \
    && git clone https://github.com/eamiddlebrook/OrthoPhyl.git \
    && cd OrthoPhyl \
    && chmod a+rw OrthoPhyl.sh \
    && cd /opt \
    && apt-get clean

ENV PATH="/opt/Path_to_gits/OrthoPhyl/:/opt/miniconda3/bin/:$PATH"
ENV Path_to_gits="/opt/Path_to_gits"
ENV OrthoPhyl="bash /opt/Path_to_gits/OrthoPhyl/OrthoPhyl.sh"
ENV ASTRAL_cmd="$Path_to_gits/ASTRAL/Astral/astral.5.7.8.jar"
ENV catfasta2phyml_cmd="$Path_to_gits/catfasta2phyml/catfasta2phyml.pl"
ENV Alignment_Assessment="$Path_to_gits/Alignment_Assessment/Alignment_Assessment_v2.py"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN echo "Container was created $(date)" > /README.txt

ENTRYPOINT ["/bin/bash", "/opt/Path_to_gits/OrthoPhyl/Ortho
