# Note
This version is copied from https://gitlab.com/sevyharris/MonteCoffee, which is a fork of https://gitlab.com/ChemPhysChalmers/MonteCoffee

I used the following commands to rewrite the git history without the 100MB+ files that were preventing automatic transfer to Github
```
git filter-branch --tree-filter 'rm -f NeighborKMC/tutorials/A_ads/result/cov_100x10.txt' HEAD
git filter-branch -f --tree-filter 'rm -f NeighborKMC/tutorials/B2_ads/result/cov_10x100.txt' HEAD
```

# What is MonteCoffee?
MonteCoffee is an application for performing kinetic Monte Carlo 
(kMC) simulations. The program is written as a set of Python modules, which
are meant to be used as a programmable application. In this sense, the user
modifies and develops on the code to customize as needed.

# Description and how to cite?
If you find MonteCoffee useful, please make a citation to:

[M. Jørgensen and H. Grönbeck, *J. Chem. Phys.* **(2018)**, *149*, 114101](https://doi.org/10.1063/1.5046635)


# Scientific articles using MonteCoffee
The following papers have used MonteCoffee:

[M. Jørgensen and H. Grönbeck, *J. Chem. Phys.* **(2018)**, *149*, 114101](https://doi.org/10.1063/1.5046635)

[M. Jørgensen and H. Grönbeck, *ACS Catal.* **(2017)**, *7*, 5054](https://doi.org/10.1021/acscatal.7b01194)

[M. Jørgensen and H. Grönbeck, *Angew. Chem. Int. Ed.* **(2018)**, *57*, 5086](https://doi.org/10.1002/anie.201802113)

[T. N. Pingel, M. Jørgensen, A. B. Yankovich, H. Grönbeck, and E. Olsson, *Nat. Commun.* **(2018)**, *9*, 2722](https://doi.org/10.1038/s41467-018-05055-1)

[M. Jørgensen and H. Grönbeck, *Top. Catal.* **(2019)**](https://doi.org/10.1007/s11244-019-01145-6)

[M. Jørgensen and H. Grönbeck, *J. Am. Chem. Soc.* **(2019)**](https://doi.org/10.1021/jacs.9b02132)

# Running a MonteCoffee simulation
The *quick_example.py* in the documentation describes a minimal working example on how to run a simulation.

for a quick-start example of CO oxidation over a nanoparticle, you can run *python test.py* 

## Documentation
Is hosted on [montecoffee.readthedocs.io](https://montecoffee.readthedocs.io/)
and is built automatically after each push.

To build documentation manually, see *documentation/HOWTO_BUILD_DOCUMENTATION.txt*

## Developers

- Mikkel Jørgensen, Chalmers University of Technology, Sweden.

- Noemi Bosio, Chalmers University of Technology, Sweden.

- Elisabeth M. Dietze, Chalmers University of Technology, Sweden.

## Troubles / Ideas / Questions?

Please post an issue to the git, and we will try to reply ASAP.
