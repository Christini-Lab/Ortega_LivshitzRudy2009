# Livshitz Rudy 2009 Model - Dynamic Clamp Optimized
---

### Introduction
This version of the Livshitz Rudy 2009 guinea pig ventricular model has been
modified for real-time performance in a dynamic clamp setting.

This repository is intendent to be used as a subtree. Follow the directions
below to add/update the model to your repository.

### Subtree installation:

Enter your local Git project
```sh
cd <PROJECT DIRECTORY>
```

Add remote URL of the model to your local project.
  * **-f** - fetch from remote immediately
  * **LivR2009_remote** - Name of remote, you can change this

```sh
git remote add -f LivR2009_remote https://pbtech-vc.med.cornell.edu/git/fro2002/dynclamp_livr2009_model.git
```

Add genetic algorithm as a subtree of the project.
  * **--prefix LivRudy2009/** - prefix denotes the directory you wish to
  put the model in, you can change this
  * **LivR2009_remote** - Model remote repository set earlier
  * **--squash** - merges all commits into one for cleaner history

```sh
git subtree add --prefix LivRudy2009/ LivR2009_remote master --squash
```

### Updating the Genetic_Algorithm subtree:
Fetch any new changes, then pull changes into subtree directory.
  * **--prefix LivRudy2009/** - Specify the model directory
  * **LivR2009_remote** - Model remote repository you will pull changes from
  * **master** - Pulling in master branch changes, you can specify another
  branch if required

```sh
git subtree pull --prefix LivRudy2009/ LivR2009_remote master --squash
```