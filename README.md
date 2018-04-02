# Livshitz Rudy 2009 Model - Ortega Modified Version
---

### Introduction
This version of the Livshitz Rudy 2009 guinea pig ventricular model. The model
is split into several branches for specific purposes. For example, the real-time
optimized version used in the 2017 Devenyi and Ortega (cited below) is located
in the RealTime branch. Underlying formulations are the same as original 2009
paper (cited below). Modifications surround formatting and implementation.

Branches:
  * **Adaptive** - Uses an adaptive time step for integration.
  * **Optimized** - Most recent real-time optimized version.
  * **RealTime** - Version used in the 2017 Devenyi and Ortega paper.
  * **Timed** - Fixed time step.
  * **master** - Current used version, similar to Optimized.

Devenyi, Ryan A., et al. "Differential roles of two delayed rectifier potassium currents in regulation of ventricular action potential duration and arrhythmia susceptibility." The Journal of physiology 595.7 (2017): 2301-2317.
Livshitz, Leonid, and Yoram Rudy. "Uniqueness and stability of action potential models during rest, pacing, and conduction using problem-solving environment." Biophysical journal 97.5 (2009): 1265-1276.

This repository can be used as a subtree. Follow the directions below to
add/update the model to your repository.

### Running Simulations
Simulation protocols are located in the 'exe' directory. Compile and run the
corresponding executable to start simulation.

### Subtree installation:

Enter your local Git project
```sh
cd <PROJECT DIRECTORY>
```

Add remote URL of the model to your local project.
  * **-f** - fetch from remote immediately
  * **LivR2009_remote** - Name of remote, you can change this

```sh
git remote add -f LivR2009_remote https://github.com/francis-ortega/Ortega_LivshitzRudy2009.git
```

Add model as a subtree of the project.
  * **--prefix LivRudy2009/** - prefix denotes the directory you wish to
  put the model in, you can change this
  * **LivR2009_remote** - Model remote repository set earlier
  * **--squash** - merges all commits into one for cleaner history

```sh
git subtree add --prefix LivRudy2009/ LivR2009_remote master --squash
```

### Updating the model subtree:
Fetch any new changes, then pull changes into subtree directory.
  * **--prefix LivRudy2009/** - Specify the model directory
  * **LivR2009_remote** - Model remote repository you will pull changes from
  * **master** - Pulling in master branch changes, you can specify another
  branch if required

```sh
git subtree pull --prefix LivRudy2009/ LivR2009_remote master --squash
```