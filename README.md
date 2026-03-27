# Exercises_and_project_work_WW

This repository contains my solutions to course exercises as well as my project work. It serves as a structured record of my learning progress and practical implementations throughout the course.

---

## Repository Structure

```text
.
├── exercises/
│   ├── day_01/
│   ├── day_02/
│   ├── day_03/
│   ├── day_04/
│   └── day_05/
│
├── project/
│   ├── 01_old_code/
│   └── wiebke_qpcr_project/
│
└── README.md
```

---

## Exercises

The `exercises/` folder contains solutions to daily assignments.
Each day has its own subfolder (e.g., `day_01`, `day_02`, etc.), allowing for a clear and chronological organization of tasks and progress.

---

## Project

The `project/` folder contains the main project developed during the course.

### qPCR Analysis Pipeline

The primary project is located in:

```
project/wiebke_qpcr_project/
```

This project implements a **modular and reproducible RT-qPCR analysis pipeline**, including:

* preprocessing and quality control
* technical replicate summarization
* normalization using reference genes
* optional inter-plate calibration
* statistical testing
* automated plotting and reporting

The goal is to move from one-off analysis scripts toward a **flexible, scalable, and reusable analysis workflow** suitable for real research applications.

This project builds on code originally developed during my master’s work and extends it into a more structured and robust pipeline.

See the project-specific README for full details:

```
project/wiebke_qpcr_project/README.md
```

---

## Purpose

* Practice and deepen understanding of Python programming
* Maintain a structured overview of completed exercises
* Apply concepts in a larger, real-world scientific project
* Develop a reusable and professional data analysis pipeline
* Document learning progress over time

---

## Getting Started

To use this repository locally:

1. Clone the repository

```bash
git clone <your-repository-url>
```

2. Navigate into the directory

```bash
cd Exercises_and_project_work_WW
```

3. Explore the contents

* `exercises/` → daily tasks
* `project/` → main project work

---

## Requirements

This repository does not enforce a single global environment.

Each project or exercise may define its own requirements.

For the qPCR pipeline, see:

```
project/wiebke_qpcr_project/environment.yml
```

or

```
project/wiebke_qpcr_project/requirements.txt
```

---

## Notes

* This repository reflects a learning process and includes exploratory as well as structured code
* Different approaches are used as new concepts are introduced throughout the course
* The project component focuses on applying programming skills to a real scientific problem

---

## Author

**Wiebke Weber**
PhD Candidate, Uppsala University




