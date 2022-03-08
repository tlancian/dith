## Discovering Polarization Niches via Dense Subgraphs with Attractors and Repulsers

[Adriano Fazzone](https://scholar.google.it/citations?user=ivW-SnEAAAAJ&hl=en) (Sapienza University, Rome), [Tommaso Lanciano](https://phd.uniroma1.it/web/LANCIANO-TOMMASO_nP1661409_EN.aspx) (Sapienza University, Rome), [Riccardo Denni](https://phd.uniroma1.it/web/RICCARDO-DENNI_nP1893279_EN.aspx) (Sapienza University, Rome), [Charalampos Tsourakakis](https://tsourakakis.com/) (ISI Foundation, Turin and Boston University, Boston), [Francesco Bonchi](http://www.francescobonchi.com/) (ISI Foundation, Turin and Eurecat, Barcelona)

<p align="center">
  <img width="600" height="300" src="https://github.com/tlancian/dith/blob/main/figure1.png">
</p>

---

_Detecting niches of polarization in social media is a first step towards deploying mitigation strategies and avoiding radicalization. In this paper, we model polarization niches as close-knit dense communities of users, which are under the influence of some well-known sources of misinformation,  and isolated from authoritative information sources. Based on this intuition we define the problem of finding a subgraph that maximizes a combination of density, proximity to a small set of nodes A (Attractors), and distance from another small set of nodes R (Repulsers). Deviating from the literature, we do not exploit text mining or sentiment analysis, nor we track the propagation of information: we only exploit the network structure and the background knowledge about the sets A and R, which are given as input. We build on very recent algorithmic advances in supermodular maximization \[Chekuri et al., SODA 22\] to provide an iterative greedy algorithm, dubbed Down in the Hollow (DITH), that converges fast to a near optimal solution. Thanks to a theoretical upper bound, we equip DITH with a practical device that allows to terminate as soon as a solution with a user-specified approximation factor is found, making our algorithm very efficient in practice.  Our experiments on very large networks confirm that our algorithm always returns a solution with an approximation factor better or equal to the one specified by the user, and it is scalable. Our case-studies in polarized settings, confirm the usefulness of our algorithmic primitive in detecting polarization niches._

---

### Requirements

The code has been tested with Java 11.


### Folders

* **datasets**: 
* **sw**: 
* **output**: 

### Execution


#### Positional arguments:



#### Optional arguments:

  	
#### Examples:



### Contacts
Mail to [fazzone@diag.uniroma1.it](mailto:fazzone@diag.uniroma1.it) and [tommaso.lanciano@uniroma1.it](mailto:tommaso.lanciano@uniroma1.it) for any question.
