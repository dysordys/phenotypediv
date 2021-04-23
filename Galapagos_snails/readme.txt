Description of snaildata.csv
============================

Rows are individual shells, so the dataset includes 1787 individuals.

Columns:

- species: Most of these are published species names, although a few are currently unpublished (e.g., "corini", "krameri MS", and "phanaticus MS"), or are morphospecies that we think are good biological species. Lineages corresponding to the type specimens are named normally (e.g., "albemarlensis"), while secondary lineages are named with an extra label (like "albermarlensis CA"). In this case, the label "CA" shows that the specimens are from a different volcano - Cerro Azul, in this case - than the type specimen. Another example is "perrus1", a population sampled from one side of Fernandina Island, while "perrus" (without any labels) is from the opposite side and represents a separate lineage.

- island: The island each individual comes from - or, in the case of the Isabela individuals, the particular volcano (AL: Alcedo Volcano, Isabela; CA: Cerro Azul, Isabela; CH: Champion islet, near Floreana; DA: Darwin Volcano, Isabela; ED: Eden islet, near Santa Cruz; ES: Espanola; FA: Fernandina; FL: Floreana; GA: Gardner islet, near Floreana; PI: Pinzon; SA: Santiago; SC: Santa Cruz; SF: Santa Fe; SL: San Cristobal; SN: Sierra Negra, Isabela; WF: Volcan Wolf, Isabela).

- site: This is either the station within which we had found the specimen (e.g. SL05-10), or is the museum lot that the shell came from (e.g. MCZ 278892).

- vegzoneHumid/vegzoneTransition/vegzoneArid: These columns correspond to the vegetative zones each species tends to be found within. For example, achatellinus is typically a humid zone species, while calvus can be found in both the transition zone and the arid zone.

- microhabitat: Each species can be generally classified as arboreal (found on vegetation off the ground) or terrestrial (found under cover or in the leaf litter on the ground), although there is often broad overlap when they are especially active (e.g., after a rain).

- Centroidsize: This measure of shell size is quantified as the average distance between each landmark used for geometric morphometrics to the centroid of those landmarks.

- shape1-shape28: These columns together characterize the locations of each individual in shell shape space. This is quantified by taking all shell landmark coordinates; aligning them such that variation due to orientation, position, and scale have been stripped away; and projected into tangent space.
