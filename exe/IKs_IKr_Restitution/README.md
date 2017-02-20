# Restitution Portrait Simulation - IKs/IKr Ratio
---

### Introduction
This protocol is designed to test the impact of IKs/IKr ratio.

Performs a restitution portrait simulation using the Livshitz Rudy 2009 model.
Relevant parameters are the start and ending BCL desired, number of beats per
BCL, and the number of calculated APDs saved. If the model crashes or the number
of stimulations does not match the number of APs detected, the protocol will
end.

### Scripts

  * **IKs_IKr_Restitution** - Script to run all ratios on a local machine. This
  script is intended to be an example, and is not realistic to run due to its
  long computation time.
  * **IKs_IKr_Restitution_Cluster_Array** - Script for the PBTech compute
  cluster. This is intended to be run as an array job, where a single core
  runs a simulation on a single ratio.