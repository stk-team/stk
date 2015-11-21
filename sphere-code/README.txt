====Disclaimer====
THIS SPHERICAL ANALYSIS CODE BASE IS CREATED AND MAINTAINED BY GURPRIT SINGH. 
ALL RIGHTS RESERVED. Copyright 2014-15.
All rights are given for usage of the code for research purposes. Users can add/modify 
the code according to their comfort. However, while publishing please acknowledge 
the author(s).
==================


Where are the files with main() function ?
- All .cxx files are in demos/<> folder.

How can I execute a demo ?
- To execute a demo:
 --- Open CMakeLists.txt and go to the end of the file
 --- At the end of the file, write the name of the demo.cxx file you want to execute (without .cxx extension).
 --- Now close CMakeLists.txt file and run following commands:
 --- cmake .
 --- make

Where are the results ?
- All results are in the folder FOLDER_WITH_CMAKELISTS.TXT_FILE/results/
- results are stored automatically with the date.
- All results are automatically categorized in datafiles, images, graphs section.
- You can use these prefixes while giving name to your file you plan to write.

