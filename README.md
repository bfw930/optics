
# UV SIR LIBS
> Sonic Information Retreival Libraries for Uncanny Valley


Summary:

    uvsir analysis libraries, purely functional



#### Audio Feature Graph Database
---
__Audio Processing and Feature Extraction__


Structure:

    core: orchestration functions for building database, importing audio data, processing and extracting features, and
    performing relational analysis, as well as data extraction (audio segments for playback)

    database: graph database for storage of source data, calculation results, segmentation pointers, extracted features,
    with search by node class or properties

    process: raw audio processing for basic transformation (constant-q transform, chromagram, etc.), temporal similarity
    segmentation, feature calculation

    analysis: feature selection and aggregation, dimensionality reduction (UMAP), similarity clustering (HDBSCAN*), and
    training of learning algorithms, as well as relationship investigation helpers with visualisation



[//]: # (this is a comment)



## Installation

Depends:

    See requirements.txt for complete list of python package dependencies

    Recommend using a python virtual environment within package root directory

    Notebook templates reference relative path to data folder within package root directory, update as required


## Usage

Usage details:

    See jupyter notebooks in notebooks folder for operational templates, usage documentation therein

    primary deps:

        numpy==1.16.4
        pandas==0.24.2

        scipy==1.3.0
        librosa==0.6.3

        scikit-learn==0.21.2
        hdbscan==0.8.22
        umap-learn==0.3.9

        matplotlib==3.1.0
        bokeh==1.2.0

