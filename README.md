## Discovering Polarization Niches via Dense Subgraphs with Attractors and Repulsers

[Adriano Fazzone](https://scholar.google.it/citations?user=ivW-SnEAAAAJ&hl=en) (Sapienza University, Rome), [Tommaso Lanciano](https://phd.uniroma1.it/web/LANCIANO-TOMMASO_nP1661409_EN.aspx) (Sapienza University, Rome), [Riccardo Denni](https://phd.uniroma1.it/web/RICCARDO-DENNI_nP1893279_EN.aspx) (Sapienza University, Rome), [Charalampos Tsourakakis](https://tsourakakis.com/) (ISI Foundation, Turin and Boston University, Boston), [Francesco Bonchi](http://www.francescobonchi.com/) (ISI Foundation, Turin and Eurecat, Barcelona)

<p align="center">
  <img width="600" height="300" src="https://github.com/tlancian/dith/blob/main/figure1.png">
</p>

---

_Detecting niches of polarization in social media is a first step towards deploying mitigation strategies and avoiding radicalization. In this paper, we model polarization niches as close-knit dense communities of users, which are under the influence of some well-known sources of misinformation,  and isolated from authoritative information sources. Based on this intuition we define the problem of finding a subgraph that maximizes a combination of density, proximity to a small set of nodes A (Attractors), and distance from another small set of nodes R (Repulsers). Deviating from the literature, we do not exploit text mining or sentiment analysis, nor we track the propagation of information: we only exploit the network structure and the background knowledge about the sets A and R, which are given as input. We build on very recent algorithmic advances in supermodular maximization \[Chekuri et al., SODA 22\] to provide an iterative greedy algorithm, dubbed Down in the Hollow (DITH), that converges fast to a near optimal solution. Thanks to a theoretical upper bound, we equip DITH with a practical device that allows to terminate as soon as a solution with a user-specified approximation factor is found, making our algorithm very efficient in practice.  Our experiments on very large networks confirm that our algorithm always returns a solution with an approximation factor better or equal to the one specified by the user, and it is scalable. Our case-studies in polarized settings, confirm the usefulness of our algorithmic primitive in detecting polarization niches._

---

### Requirements

The code has been implemented in Java 11. The code runs Algorithm 2, assuming that proximities and distances are pre-computed and given in input.


### Folders

* **datasets**: folder containing the graph and an example of the pre-computed proximity/distance values for a specific instance of Problem 1. Due to Github single-file space limits, we provide here only a sample of the Datasets. All the graphs employed in this work can be found [here](https://www.dropbox.com/sh/03ouxvlv4jjvht9/AAClYn58iD3wCIrGjdIxLbjDa?dl=0).
* **sw**: folder containing all the software necessary to run DITH.
* **output**: folder that contains the output files.

### Execution
To run the code go in ./sw/bin and run the following command:
java -Xms4g -Xmx4g algorithms.SuperGreedyPP algorithm data_directory graph_file_name set_A set_B lambda_1 lambda_2 complete_proximity_file_name complete_distance_file_name T gamma output_data_directory experiment_id 

#### Positional arguments:

* **algorithm**: algorithm to be executed. Possible choices are {"DITH", "DITH-1"}.
* **data_directory**: directory containing the input graph.
* **graph_file_name**: name of the file representing the input graph.
* **set_A_id**: string representation of set A (ignored by the software and reported only in output).
* **set_B_id**: string representation of set B (ignored by the software and reported only in output).
* **lambda_1**: lambda_1 coefficient.
* **lambda_2**: lambda_2 coefficient.
* **complete_proximity_file_name**: complete file name of the file containing all the proximity values.
* **complete_distance_file_name**: complete file name of the file containing all the distance values.
* **T**: maximum number of iterations the algorithm is allowed to do.
* **gamma**: desired approximation factor.
* **output_data_directory**: directory that will contains the output.
* **experiment_id**: string identifier that represents the performed experiment.

  	
#### Examples:

java -Xms4g -Xmx4g algorithms.SuperGreedyPP "DITH" ../../datasets/gunsense gunsense.tsv "{'1548714900', '2551313612', 'fabaceae', 'gregmercer1', 'pppatticake'}" "{'2984727962', 'koopac7', 'rickcanton'}" 196.23 242.95 ../../datasets/gunsense/gunsense_proximity.tsv ../../datasets/gunsense/gunsense_distance.tsv 10000 "0.99" "../../output" "gunsense_TEST_experiment" 

java -Xms4g -Xmx4g algorithms.SuperGreedyPP "DITH" ../../datasets/vaxnovax vaxnovax.tsv "{'RobertoBurioni', 'nzingaretti'}" "{'DiegoFusaro', 'GiorgiaMeloni'}" 51.7 56.4 ../../datasets/vaxnovax/vaxnovax_proximity.tsv ../../datasets/vaxnovax/vaxnovax_distance.tsv 10000 "0.99" "../../output" "vaxnovax_TEST_experiment"


### Contacts
Mail to [fazzone@diag.uniroma1.it](mailto:fazzone@diag.uniroma1.it) and [tommaso.lanciano@uniroma1.it](mailto:tommaso.lanciano@uniroma1.it) for any question.
