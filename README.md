# local_descriptors_object_recognition

This repository contains an object recognizer based on local descriptors, developed by the year 2005.

It implements those papers:

`Loncomilla, P., Ruiz-del-Solar, J. (2005). Improving SIFT-Based Object Recognition for Robot Applications. In: Roli, F., Vitulano, S. (eds) Image Analysis and Processing – ICIAP 2005. ICIAP 2005. Lecture Notes in Computer Science, vol 3617. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11553595_133`

`Loncomilla, P., Ruiz-del-Solar, J. (2006). A Fast Probabilistic Model for Hypothesis Rejection in SIFT-Based Object Recognition. In: Martínez-Trinidad, J.F., Carrasco Ochoa, J.A., Kittler, J. (eds) Progress in Pattern Recognition, Image Analysis and Applications. CIARP 2006. Lecture Notes in Computer Science, vol 4225. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11892755_72`

`Loncomilla, P., Ruiz-del-Solar, J., Martínez, L. (2016). Object recognition using local invariant features for robotic applications: A survey. Pattern Recognition,Volume 60, Pages 499-514, ISSN 0031-3203, https://doi.org/10.1016/j.patcog.2016.05.021`


This detector rejects wrong detection hypothesis by applying a cascade of several tests.

Note: This code was written using a very old C++ standard / compiler. However, it was slighty modified for compiling in current compilers, but it generates a lot of warnings.

An example use is shown in main.cpp


Example inputs:

Reference image:

![Reference image](https://github.com/ploncomi/local_descriptors_object_recognition/blob/main/uch010a.jpg)

Test image:

![Test image](https://github.com/ploncomi/local_descriptors_object_recognition/blob/main/uch010b.jpg)

Result:

![Result](https://github.com/ploncomi/local_descriptors_object_recognition/blob/main/imdraw.jpg)
