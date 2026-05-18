# Truss Analysis and Deflection Application Using Python

## CE257 Capstone Project – 2026

This project was developed as part of the CE257 Computer Programming course in the Department of Civil Engineering at Kwame Nkrumah University of Science and Technology (KNUST).

The application performs structural analysis of two-dimensional statically determinate trusses using the Method of Joints and computes truss deflections using the Unit Load Method.

---

## Project Overview

The software allows users to:

- Draw truss systems interactively
- Create joints and members
- Apply external loads
- Define pin and roller supports
- Analyze truss member forces
- Compute support reactions
- Calculate truss deflections
- Display the deformed shape graphically

The application includes an interactive graphical user interface (GUI) built using Tkinter and graphical visualization using Matplotlib.

---

## Features

- Interactive truss drawing interface
- Adjustable grid system
- Member connectivity creation
- Load application at joints
- Pin and roller support assignment
- Method of Joints analysis
- Force visualization using color coding
- Deflection analysis using the Unit Load Method
- Numerical and graphical result display
- Deformed shape visualization

---

## Structural Analysis Methods

### Method of Joints

The application analyzes statically determinate trusses using equilibrium equations:

- SFx = 0
- SFy = 0

The solver computes:
- Internal member forces
- Support reactions

---

### Unit Load Method

The Unit Load Method is used to determine:
- Horizontal joint deflections
- Vertical joint deflections
- Deformed truss shape

---

## Technologies and Libraries Used

- Python 3
- Tkinter
- NumPy
- Matplotlib
- Math Library

---

## Requirements

To run this application, ensure the following are installed:

- Python 3.x
- NumPy
- Matplotlib
- Tkinter

The application can be executed using any Python IDE or terminal environment.

---

## Installation

Install the required libraries using:

```bash
pip install numpy matplotlib
```

---

## Running the Application

Open the project folder in any Python IDE or terminal and run:

```bash
python main.py
```

---

## Project Structure

```text
project-folder/
¦
+-- main.py
+-- README.md
+-- report.pdf
+-- screenshots/
```

---

## User Guide

1. Launch the application
2. Set grid spacing
3. Add joints
4. Connect joints using members
5. Apply loads
6. Define support conditions
7. Click the SOLVE button
8. View:
   - Member forces
   - Support reactions
   - Deformed shape

---

## Verification and Testing

The application was tested using several statically determinate truss configurations.

Manual calculations using the Method of Joints were compared with the software results, and complete agreement was obtained, confirming the accuracy and reliability of the program.

---

## Group Members

| Name | Index Number |
|---|---|
| Elvis Aseye Fasemkye | 3948824 |
| Blessed Isaac Andoh | 3947724 |
| Ishmael Oduro | 3950324 |

---

## Department

Department of Civil Engineering  
Kwame Nkrumah University of Science and Technology (KNUST)

---

## Academic Year

May 2026

---

## Future Improvements

Possible future enhancements include:

- Analysis of statically indeterminate trusses
- Stress and strain calculations
- File saving and loading functionality
- Improved graphical interface
- Animated deformation visualization

---

## License

This project was developed strictly for academic and educational purposes.

---

## Repository

Hosted on GitHub.