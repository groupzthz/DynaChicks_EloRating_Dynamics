# DynaChicks_EloRating_Dynamics
### Chicken hierarchy evaluations of large and small groups across maturation

Data and analysis scripts for the paper: "Social structure and interactions of laying hens..."

## Data files:

### used for **Analyses_Interactions_Structure.R**:

1. **Individuals.csv**  
   Contains all existing individuals in alphabetical order for each pen  

   | **Column** | **Description** |
   |------------|-----------------|
   | **Number** | Counter from 1 to N (N is the group size of the pen) |
   | **Pen**    | Animal group from A to F (6 groups) |
   | **ID**     | Identifier of an individual within a group (only unique within group, not between groups) |

2. **MatureInteractionsFull.csv**  
   Contains all observed social interactions at 24 WoA   

   | **Column**     | **Description** |
   |----------------|-----------------|
   | **Date**       | Observation date |
   | **Pen**        | Animal group from A to F (6 groups) |
   | **Winner**     | Identifier of the winner of the interaction |
   | **Loser**      | Identifier of the loser of the interaction |
   | **Condition**  | Situation in which the animals were observed (Normal, HQ or Feeder) |
   | **VideoTime**  | Time point within video recording (not related to the actual time) |
   | **Code**       | Type of agonistic interaction observed (Avoidance, Peck, Threat, or Fight) |
   | **VideoLocation** | Indicator for large groups from which camera the video came |

### used for **Fear_Recognition_Badges.R**:

3. **ArenaTest\*.csv**  
   Contains the data collected for the fear (Test 1) and the recognition (Test 2) test, as well as other parameters for the birds.  
   Only columns used in the present study are listed. If they exist in only one of the two files, they are numbered.   

   | **Column**          | **Description** |
   |---------------------|-----------------|
   | **Pen**             | Animal group from A to F (6 groups) |
   | **ID**              | Identifier of an individual within a group (only unique within group, not between groups) |
   | **Side/Side_Flockmate** | Location of the conspecific with respect to the starting position (A/S or Right/Left) |
   | **FirstChoice**     | Whether the hen first entered the field in proximity to the conspecific or the feed/the flockmate or the non-flockmate, left empty or NA if the hen made no choice (Social/Feed/ " " or Flockmate/NonFlockmate/NA) |
   | **Weight/Weight26** | Body mass of the animal (g) |
   | **Box (1)**         | Latency to emerge from the box (s) |
   | **Feeding (1)**     | Boolean indicating whether the hen ate (0/1) |
   | **FirstEntry (2)**  | In rare occasions, the chicken escaped when placed in the arena and entered a vicinity during fleeing from humans; those were not counted as choices but as first entries and excluded from the choice-based analysis. |
   | **Flockmate (2)**   | Total time spent in vicinity of the flockmate |
   | **Non-Flockmate (2)** | Total time spent in vicinity of the non-flockmate |

4. **ThermalResults.csv**  
   Contains the comb size information needed for the badges of status analysis. Other data in the sheet is not used.  

   | **Column**    | **Description** |
   |---------------|-----------------|
   | **Pen**       | Animal group from A to F (6 groups) |
   | **ID**        | Identifier of an individual within a group (only unique within group, not between groups) |
   | **CombSize1** | Comb size measured at 16 WoA (in cm²) (not used for analysis) |
   | **CombSize2** | Comb size measured at 26 WoA (in cm²) |

5. **IndividualsOutput.csv**  
   File created within "Analyses_Interactions_Structure.R".  
   Contains an overview of the data collected on an individual level.  

   | **Column**        | **Description** |
   |-------------------|-----------------|
   | **ID**            | Identifier of an individual within a group (only unique within group, not between groups) |
   | **Wins**          | Sum of all interactions the animal won |
   | **Losses**        | Sum of all interactions the animal lost |
   | **elos**          | Randomized Elo rating score of the animal |
   | **rank**          | Rank of the animal based on Elo rating (1 = highest rank) |
   | **sum**           | Total number of interactions an animal was involved in |
   | **scaleElos**     | Scaled Elo rating (scaled within pen) |
   | **physAggr**      | Amount of high-intensity aggression the animal dealt |
   | **nonphysAggr**   | Amount of low-intensity aggression the animal dealt |
   | **HQ**            | Amount of aggression the animal dealt during HQ |
   | **Feed**          | Amount of aggression the animal dealt during FC |
   | **Normal**        | Amount of aggression the animal dealt without resources present (Control) |
   | **physAggrRec**   | Amount of high-intensity aggression the animal received |
   | **nonphysAggrRec** | Amount of low-intensity aggression the animal received |
   | **HQRec**         | Amount of aggression the animal received during HQ |
   | **FeedRec**       | Amount of aggression the animal received during FC |
   | **NormalRec**     | Amount of aggression the animal received without resources present (Control) |
   | **Pen**           | Animal group from A to F (6 groups) |
   | **Group_Size**    | Group size (small/large) |

## R scripts:

1. **Analyses_Interactions_Structure.R**  
   Calculation of the Elo ratings for the chickens and analysis of hierarchy characteristics, interactions, and other factors to compare small and large group sizes.  
   Creates the data file `IndividualsOutput.csv`, used in **Fear_Recognition_Badges.R**.

2. **helper_functions.R**  
   Contains functions for the Elo rating calculations and other computations.

3. **Fear_Recognition_Badges.R**  
   Analysis of badges of status, the fear test, and the recognition test.