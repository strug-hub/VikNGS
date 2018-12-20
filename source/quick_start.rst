.. _quick_start:

Quick Start
==================================

Installing Software
------------------------------

To download VikNGS, go to the `GitHub repository <https://github.com/ScottMastro/vikNGS/tree/master/>`_ and press "Clone or download" and Download the ZIP version of the repository as shown below.

.. figure:: resources/download.png

User Interface
------------------------------

The Release folder of the `VikNGS GitHub repository <https://github.com/ScottMastro/vikNGS/tree/master/>`_ contains precompiled versions of VikNGS for Windows, Mac and Linux operating systems.

Running the file VikNGS-X.X.X in the appropriate operating system folder should start up the user interface. The available versions were compiled on the following systems, so please try to run it under a similar setting:

* Windows: Windows 10 64x (compiled with Microsoft Visual C++ 2017)

* Mac: macOS 10.13.4 High Sierra

* Linux: Ubuntu 18.04 64x

If there is an issue running the software, we recommend trying a different system or  :ref:`compiling the software from the source code<source_code>`.

Running Example on User Interface
------------------------------

In the same directory where the VikNGS-X.X.X application lives, two example files are present. On the Association Test tab of VikNGS, provide *example.vcf* as input in the section labelled "VCF File" and *example_info.txt* as input in the section labelled "Sample Information File". Clicking the "RUN" button should then trigger the association tests using these files as input.


Running Command Line
------------------------------

From the command line, run the following commands:

| ``wget -O vikNGS https://github.com/ScottMastro/vikNGS/archive/master.zip``
| ``unzip vikNGS``
| ``cd vikNGS/bin``
| ``make``

