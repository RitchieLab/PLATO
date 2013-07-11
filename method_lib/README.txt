The files in this folder should be organized as follows:

-root (method_lib):
Should contain the Makefile.am and the following abstract base classes:
  - Method
  - AnalysisMethod
  - FilterMethod
  - InputMethod
  - OutputMethod

- util:
Should contain any classes that will be widely used (like Sample, Marker, etc)

- analysis:
Should contain concrete subclasses of AnalysisMethod
These classes perform some sort of analysis on the data

- filters:
Should contain concrete subclasses of FilterMethod
These classes filter the data

- input:
Should contain concrete subclasses of InputMethod
These classes take input data

- output:
Should contain concrete subclasses of OutputMethod
These classes output (reformat) data
