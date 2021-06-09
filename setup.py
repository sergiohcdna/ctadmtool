from setuptools import setup , find_packages
import pathlib

here = pathlib.Path( __file__ ).parent.resolve()

long_description = ( here / 'README.md' ).read_text( encoding='utf-8' )

setup(

    name='ctaAnalysis' ,

    version='0.1.1' ,

    description='A tool for indirect dark matter searches with CTA' ,

    long_description=long_description ,

    long_description_content_type='text/markdown' ,

    author='Sergio Hernández Cadena, Judit Pérez Romero, Miguel Sánchez Conde' ,

    author_email='skerzot@ciencias.unam.mx' ,

    classifiers=[
        'Development Status :: 4 - Beta' ,
        'Intended Audience :: Science/Research' ,
        'License :: OSI Approved :: MIT License' ,
        'Programming Language :: Python :: 3 :: Only' ,
        'Programming Language :: Python :: 3.6'
        ] ,

    keywords=[
        'cta' , 'dark matter' , 'indirect searches' , 'science' ,
        'indirect dm searches'
        ] ,

    #   packages=find_packages( where='ctadmtool' )
    packages=find_packages( exclude=( 'old' , ) ) ,

    install_requires=[ 'scipy' , 'numpy' , 'matplotlib' , 'ebltable', 'tqdm' ] ,

    include_package_data=True,

    package_data={ '' : [ 'data/*.dat' , 'pfiles/*.par' , 'pfiles/*.txt'  ] }

)

