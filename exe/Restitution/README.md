# Restitution Portrait Simulation
---

### Introduction
Performs a restitution portrait simulation using the Livshitz Rudy 2009 model.
Relevant parameters are the start and ending BCL desired, number of beats per
BCL, and the number of calculated APDs saved. If the model crashes or the number
of stimulations does not match the number of APs detected, the protocol will
end.

### Run Instructions

Compile executable
```sh
make
```

Run simulation
```sh
./Restitution_Portrait
```