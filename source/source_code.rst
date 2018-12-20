.. _source_code:

Source Code
==================================

Compiling the VikNGS User Interface
------------------------------

The VikNGS source code is contained in ``vikNGS/src/`` and the files specific to the graphical user interface (GUI) are found in ``vikNGS/src/gui``. To build the GUI version of the software, we recommend downloading and installing `QT 5.11+ and QT Creator <http://doc.qt.io/qt-5/index.html>`_. 

Open QT Creator after downloading and press the "Open Project" button to load the user interface QT project. Navigate to the directory where the vikNGS was downloaded and load the file ``vikNGS/src/gui/gui.pro``. This should load the source code and prompt you to choose a compiler. After selecting a compiler, the program can be build by switching to "Release" mode and pressing the top green arrow as seen below:

.. figure:: resources/qt_build.png
   :target: source_code.html
   :alt: QT build button

Compiling the VikNGS Command Line Tool
------------------------------

A Makefile is provided to compile the code for command line use. The command line-specific files are contained in ``vikNGS/src/cmd`` 

From the command line, run the following commands:

| ``wget https://github.com/ScottMastro/vikNGS/archive/master.zip``
| ``unzip master.zip``
| ``cd vikNGS-master/bin``
| ``make``

To test the binary executable file, try running the following command:

| ``./vikNGS --vcf example.vcf --sample example_info.txt``
