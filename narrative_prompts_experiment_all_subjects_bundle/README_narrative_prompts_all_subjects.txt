Narrative LLM prompts for the experimental portfolio-choice dataset (all subjects).

Each subject has 50 trials. The prompt reveals trials t=1..49 (prices + dollar spends) and asks the model to predict t=50 allocation given the real t=50 prices.

Stock names:
  X -> PapayaTech
  Y -> AxolotlWorks
  Z -> QuipuQuantum

Files:
  - narrative_prompt_subject<ID>.txt  (prediction prompt for that subject)
  - narrative_answer_subject<ID>.txt  (ground truth for t=50 for that subject)
  - narrative_index_all_subjects.csv  (one row per subject with t=50 prices and ground truth spends)