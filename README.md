# Overlapping Community Detection: HCP Structural Connectome

Code to apply and assess the results of overlapping community-detection algorithms on a right-hemisphere structural connectome estimated from Human Connectome Project data.

Depends on a toolbox for fitting benchmark overlapping networks and comparing the performance of algorithms on these networks: [here](https://github.com/NeuralSystemsAndSignals/OverlappingCommunityDetection).

Main results can be reproduced as follows:

```matlab
Fig2
```

Comments:
The scripts to reproduce figures point to the source codes (if any) that maybe used to generate any result used in the respective script. All such source codes are in `Peripheral` folder.

Other analyses not presented in the form of figures:
1. Jaccard Coefficient Analysis: Can be reproduced by running `jaccard_test` in `Peripheral`.
2. F1 and Omega score over synthetic benchmarks: Can be reproduuced using `compute_metrics` in `Peripheral`.



