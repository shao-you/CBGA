CBGA
====

Clustering–Based Genetic Algorithm

This program is written from 2012/8 in C++ language.


The main problem I target against is to find an 'optimized path' in a Wireless Sensor Network (WSN).
In short, use Genetic Algorithm (GA) for TSPN problem.
  (Mainly, find the approximately short path.)
WSN consists of many randomly deployed sensors with different sensor radii.

First, I aggregate all the sensors into some clusters which are disjoint mutually.
I proposed an innovative method to check wheather sensors are belongs to same cluster, which I think where the most variable part is. And then, a watpoint is given for each cluster. Mobile Robot (data mule) can communicate with all sensors of one cluster by visiting only one waypoint.

Second, run genetic algorithm to find a proper order to visit all cluster waypoints.
  (Cloest Position Algorithm is imposed among the fitness evaluation/crossover/mutation.)

Third,


Future work:
Open the program to the public in a client-server way. User around the world can tune the parameters in his/her needs and then easily run it by 'clik' button on the web browser. The user interface may be designed in a dynamic/interactive style.

