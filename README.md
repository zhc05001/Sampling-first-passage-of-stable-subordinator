# Exact sampling of the first passage of a stable subordinator
This R package is for exact sampling of the first passage event of a stable subordinator across a boundary.  The algorithms implemented in the package are developed in 

- Chi, Z. (2024). *Complexity of exact sampling of the first passage of a stable subordinator*. [arXiv:xxxxx](http://merlot.stat.uconn.edu/~zhc05001/)

Currently the package only samples the undershoot and jump at the first pasage across constant level 1.  However, using the package, the entire event, which include the time of the first passage as well as the undershoot and jum, across any constant level or non-constant regular boundary can be easily sampled using the package; see Algorithm 2.1 in Chi (2024).  A short R function is provided below.

To sample `n`


