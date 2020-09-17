# edge-ts

Example code for transforming parcellated fMRI BOLD activity into co-fluctuation time series ("edge time series") and obtaining the RSS time series.

The RSS time series is thresholded to retain high- and low-amplitude frames, which are used to obtain corresponding estimates of functional connectivity (FC).

We then cluster these FC patterns using modularity maximization.
